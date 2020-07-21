#!/usr/bin/env python

# Deal with the lsf logfile
import os
if not os.path.exists("loglsf"):
    os.makedirs("loglsf")

configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]


## Parsing output directory
OUTDIR = config["resources"]["outdir"]

REF_GENOMES, = glob_wildcards(GENOMEDIR + "/{ref}_ref.fa")
ALIGNERS=["ngmlr", "minimap2"]
MOVIES, = glob_wildcards(FASTQDIR + "/{movie}.fastq.gz")
SVCALLER=["sniffles", "pbsv", "svim"]


# Parameter: sample_name
sample = "sv_sample01"
if "sample_name" in config:
    sample = config['sample_name']

checkpoint sniffles:
    input:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_sniffles.vcf"
    params:
        sniffles= config["tools"]["SNIFFLES"],
        min_mq = 20,
        read_support = 5,
        min_length = 50
    envmodules:
        "intel/18.0.1"
    shell:
        """
        {params.sniffles} -m {input} -v {output} \
                -s {params.read_support} \
                -q {params.min_mq} \
                -l {params.min_length} \
                --genotype \
                --report_read_strands
        """

checkpoint pbsv_signature:
    input:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
        temp(OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_pbsv.svsig.gz")
    conda:
        "../envs/pacbio.yml"
    params:
        map_Q = 20
    shell:
        """
        pbsv discover -q {params.map_Q} {input} {output}
        
        """

rule pbsv_call_sv:
    input:
        svsig = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_pbsv.svsig.gz",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_pbsv.vcf"
    conda:
        "../envs/pacbio.yml"
    threads:
        10
    shell:
        """
        pbsv call -j {threads} {input.ref} {input.svsig} {output}

        """

rule svim_call:
    input:
        BAM = OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants.vcf"
    conda:
        "../envs/svim.yml"
    shell:
        """
        svim alignment --sample {wildcards.sample} --min_mapq 20 --min_sv_size 50 \
         {wildcards.ref}/{wildcards.aligner}/sv_calls/svim_{wildcards.sample}/ {input.BAM} {input.ref}
        """

rule svim_filter:
    input:
        OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants.vcf"
    output:
        vcf = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_svim.vcf"
    params:
        filtered = OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants_filtered.vcf.gz"
    shell:
        """
        cat {input}  \
        | sed 's/INS:NOVEL/INS/g' \
        | sed 's/DUP:INT/DUP/g' \
        | sed 's/DUP:TANDEM/DUP/g' \
        | awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>5) {{ print $0 }} }} }}' > {output.vcf}

         #cat {input}  \
         #| awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         #else {{ if($6>5) {{ print $0 }} }} }}' \
         #| bgzip > {params.filtered}

         #tabix {params.filtered}

         """

rule filter_vcf:
    input:
        VCF = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}.vcf",
        BED = GENOMEDIR + "/{ref}_chr.bed"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}_filtered.vcf"
    params:
        bcftools = config["tools"]["BCFTOOLS"]
    shell:
        """
        {params.bcftools} view -T {input.BED} {input.VCF} \
                | {params.bcftools} view -f PASS \
                | {params.bcftools} sort > {output}

        """

rule not_pass_sv:
    input:
        VCF = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}.vcf",
        BED = GENOMEDIR + "/{ref}_chr.bed"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}_not_pass.vcf"
    params:
        bcftools = config["tools"]["BCFTOOLS"]
    shell:
        """
        {params.bcftools} view -T {input.BED} {input.VCF} \
                | {params.bcftools} view -e'FILTER="PASS"' \
                | {params.bcftools} sort > {output}

        """

rule sviper:
    input:
        srs_bam = OUTDIR + f"/srs_bam/{{ref}}_{sample}.bam",
        lrs_bam = OUTDIR + f"/{{ref}}/{{aligner}}/alignment/{sample}_merge.bam",
        ref = GENOMEDIR + "/{ref}_ref.fa",
        vcf=  OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}_filtered.vcf"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/polished/{sample}_{svcaller}_polished.vcf"
    params:
        sviper = config["tools"]["SVIPER"],
        prefix= OUTDIR + "/{ref}/{aligner}/sv_calls/polished/{sample}_{svcaller}_polished"
    shell:
        """
        {params.sviper} -s {input.srs_bam} -l {input.lrs_bam} -r {input.ref} \
                -c {input.vcf} -o {params.prefix} --output-polished-bam

        """

rule merge_sv_filtered:
    input:
        expand(OUTDIR + "/{{ref}}/{aligner}/sv_calls/{sample}_{svcaller}_filtered.vcf", ref=REF_GENOMES, aligner=ALIGNERS, sample=sample, svcaller=SVCALLER)
    output:
        vcf= temp(OUTDIR + "/{ref}/sv_combined/filtered/tmp.vcf"),
        fofn = OUTDIR + "/{ref}/sv_combined/filtered/samples.fofn",
        slist= OUTDIR + "/{ref}/sv_combined/filtered/samples_name.txt"
    params:
        survivor = config["tools"]["SURVIVOR"]
    shell:
        """
        ls {input} > {output.fofn} ;
        awk -F '/|_' '{{print $16"_"$21}}' {output.fofn} > {output.slist};
        {params.survivor} merge {output.fofn} 1000 1 1 1 0 50 {output.vcf}
        """

rule merge_sv_polished:
    input:
        expand(OUTDIR + "/{{ref}}/{aligner}/sv_calls/polished/{sample}_{svcaller}_polished.vcf", ref=REF_GENOMES, aligner=ALIGNERS, sample=sample, svcaller=SVCALLER)
    output:
        vcf= temp(OUTDIR + "/{ref}/sv_combined/polished/tmp_polished.vcf"),
        fofn = OUTDIR + "/{ref}/sv_combined/polished/samples_polished.fofn",
        slist= OUTDIR + "/{ref}/sv_combined/polished/samples_name_polished.txt"
    params:
        survivor = config["tools"]["SURVIVOR"]
    shell:
        """
        ls {input} > {output.fofn} ;
        awk -F '/|_' '{{print $16"_"$22}}' {output.fofn} > {output.slist};
        {params.survivor} merge {output.fofn} 1000 1 1 1 0 50 {output.vcf}
        """


rule high_conf:
    input:
        vcf= OUTDIR + "/{ref}/sv_combined/filtered/tmp.vcf",
        slist= OUTDIR + "/{ref}/sv_combined/filtered/samples_name.txt"
    output:
        rehead_vcf= OUTDIR + "/{ref}/sv_combined/filtered/sv_merge.vcf",
        high_conf= OUTDIR + "/{ref}/sv_combined/filtered/high_conf_{ref}.vcf"
    params:
        bcftools = config["tools"]["BCFTOOLS"]
    shell:
        """
        {params.bcftools} reheader -s {input.slist} {input.vcf} -o {output.rehead_vcf} &&
        {params.bcftools} filter -i 'INFO/SUPP = "6"' -o {output.high_conf} {output.rehead_vcf}
        """

rule high_conf_polished:
    input:
        vcf= OUTDIR + "/{ref}/sv_combined/polished/tmp_polished.vcf",
        slist= OUTDIR + "/{ref}/sv_combined/polished/samples_name_polished.txt"
    output:
        rehead_vcf= OUTDIR + "/{ref}/sv_combined/polished/sv_merge.vcf",
        high_conf= OUTDIR + "/{ref}/sv_combined/polished/high_conf_{ref}.vcf"
    params:
        bcftools = config["tools"]["BCFTOOLS"]
    shell:
        """
        {params.bcftools} reheader -s {input.slist} {input.vcf} -o {output.rehead_vcf} &&
        {params.bcftools} filter -i 'INFO/SUPP = "6"' -o {output.high_conf} {output.rehead_vcf}
        """

rule merge_sv_nopass_vcf:
    input:
        expand(OUTDIR + "/{{ref}}/{aligner}/sv_calls/{sample}_{svcaller}_not_pass.vcf", ref=REF_GENOMES, aligner=ALIGNERS, sample=sample, svcaller=SVCALLER)
    output:
        vcf= OUTDIR + "/{ref}/sv_nopass/sv_merge_nopass.vcf",
        fofn = OUTDIR + "/{ref}/sv_nopass/samples.fofn",
        slist= OUTDIR + "/{ref}/sv_nopass/samples_name.txt"
    params:
        survivor = config["tools"]["SURVIVOR"]
    shell:
        """
        ls {input} > {output.fofn} ;
        awk -F '/|_' '{{print $16"_"$21}}' {output.fofn} > {output.slist};
        {params.survivor} merge {output.fofn} 1000 1 1 1 0 50 {output.vcf}
        """
