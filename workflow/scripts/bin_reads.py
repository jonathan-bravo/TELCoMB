#!/usr/bin/env python

from argparse import ArgumentParser
import gzip
from Bio import SeqIO


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outdir', required=True)
    return parser.parse_args()


def gen_read_length_clusters(reads, outdir):
    for read in reads:
        seq_len = len(str(read.seq))
        outfile = f'{outdir}/{seq_len}.rl.bins.fasta.gz'
        with gzip.open(outfile, "at") as o:
            SeqIO.write(read, o, 'fasta')


def main():
    args = parse_args()
    reads = (line for line in SeqIO.parse(open(args.infile, 'r'), 'fastq'))
    gen_read_length_clusters(reads, args.outdir)


if __name__ == "__main__":
    main()