#!/usr/bin/env bash

source /cluster/home/fangzi/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

module load gcc/4.8.5 python/3.6.4 

snakemake --jobs 100 --latency-wait 60 -rp --use-conda --use-envmodules --cluster-config cluster.json --cluster "bsub -J {cluster.jobname} -n {cluster.ncore} -W {cluster.jobtime} -o {cluster.logi} -R \"rusage[mem={cluster.memo}]\""
