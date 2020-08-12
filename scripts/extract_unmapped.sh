#!/usr/bin/env bash

module load gcc/4.8.2 samtools/1.8

#to check which reads only mapped to only one reference but not the other
DIR=/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/fang/Nanopore/sv_two_assembly/

for ref in Angus UCD
do
	for aligner in minimap2 ngmlr
	do
		samtools  view -b -f 4 ${DIR}${ref}/${aligner}/alignment/bov_168_merge.bam > ${ref}_${aligner}_unmapped.bam
	done
done
