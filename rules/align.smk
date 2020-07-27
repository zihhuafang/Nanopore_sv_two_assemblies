#!/usr/bin/env python

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
sample = "sv_sample"
if "sample_name" in config:
    sample = config['sample_name']

#align reads from individual sequencing runs to different genome assemblies using ngmlr
rule ngmlr_aln:
    input:
        fastq = FASTQDIR + "/{movie}.fastq.gz",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        OUTDIR + "/{ref}/ngmlr/alignment/aln_{movie}.bam"
    params:
        ngmlr = config["tools"]["NGMLR"],
        sample = config['sample_name']
    threads:
        20
    envmodules:
        #"gcc/4.8.2",
        #"samtools/1.8"
        "gcc/4.8.5",
        "samtools/1.10"
    shell:
        """
        zcat {input.fastq} | \
        {params.ngmlr} --presets ont -t {threads} -r {input.ref} --rg-id {wildcards.movie} --rg-sm {params.sample} | \
                samtools sort -@ {threads} -T $TMPDIR > {output}
        """

#align reads from individual sequencing runs to different genome assemblies using minimap2
rule pbmm2_aln:
    input:
        fastq = FASTQDIR + "/{movie}.fastq.gz",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        temp(OUTDIR + "/{ref}/minimap2/alignment/aln_{movie}.bam")
    threads:
        10
    params:
        minimap2 = config["tools"]["MINIMAP2"],
        sample = config['sample_name']
    envmodules:
        #"gcc/4.8.2",
        #"samtools/1.8"
        "gcc/4.8.5",
        "samtools/1.10"
    shell:
        """
        {params.minimap2} -ax map-ont --MD --eqx -L -O 5,56 -E 4,1 -B 5 \
         --secondary=no -z 400,50 -r 2k -Y \
         -R "@RG\\tID:{wildcards.movie}\\tSM:{params.sample}" \
         -t {threads} {input.ref} {input.fastq}  \
         | samtools sort -@ {threads} -T $TMPDIR > {output}
        """

#merge bam files
rule merge_bam:
    input:
        lambda wildcards: expand(OUTDIR + "/{{ref}}/{{aligner}}/alignment/aln_{movie}.bam",movie=MOVIES)
    output:
        BAM = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"),
        BAI = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam.bai")
    output:
        BAM = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"),
        BAI = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam.bai")
    threads: 10
    envmodules:
        #"gcc/4.8.2",
        #"samtools/1.8"
        "gcc/4.8.5",
        "samtools/1.10"
    shell:
        """
        samtools merge -@ {threads} - {input} | \
                samtools sort -@ {threads} -T $TMPDIR > {output.BAM} && samtools index -@ {threads} {output.BAM}

        """
