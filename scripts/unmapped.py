#! /usr/bin/env python
import pandas as pd
import pysam
import sys
from argparse import ArgumentParser


def get_names(names):
    with open(names, 'r') as infile:
        dat = pd.read_csv(infile,usecols=[0],sep='\t')
        n= dat["Read_name"].tolist()
    if '' in n:
        n.remove('')
    return n


def extract_reads(options):
    n = get_names(options.names)
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    out = open(options.out, 'w')
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                start=x.reference_start
                end=x.reference_end
                mq=x.mapping_quality
                #chrom=x.target_name
                chrom=x.reference_name
                rl=x.query_length
                out.write("%s %s %s %s %s %s\n" % (name, chrom, start, end, mq,rl))


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-n', '--names', help='list of read names to extract', required=True)
    parser.add_argument('-o', '--out', help='file name for extracted alignments', required=True)   
    options = parser.parse_args()
    extract_reads(options)

