# Nanopore_sv_two_assemblies

This is a snakemake pipeline for aligning cattle Oxford Nanopore data to two cattle genome assembiles (ARS-UCD1.2 and UOA_Angus) and calling structural variants.

# Dependencies
aligners:
ngmlr
minimap2

sv callers:
svim
pbsv
sniffles

survivor
bcftools
samtools
vep


Conda environment files containing the dependencies can be found in the envs directory. These can be automatically installed and managed by providing the --use-conda argument to snakemake: snakemake --use-conda ....
