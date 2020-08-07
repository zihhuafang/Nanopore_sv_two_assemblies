# Nanopore_sv_two_assemblies

This is a snakemake pipeline for aligning cattle Oxford Nanopore data to two cattle genome assembiles (ARS-UCD1.2 and UOA_Angus) using minimap2 and ngmlr and calling structural variants using different sv callers. The output is the high confidance sv set. 

# Dependencies
aligners:
ngmlr
minimap2

sv callers:
svim
pbsv
sniffles

others:
survivor
bcftools
samtools
vep


Conda environment files containing the some dependencies can be found in the envs directory. These can be automatically installed and managed by providing the --use-conda argument to snakemake: snakemake --use-conda ....
