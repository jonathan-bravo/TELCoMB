#!/usr/bin/env python

from argparse import ArgumentParser
from Bio import SeqIO


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--duplicates', required=True)
    parser.add_argument('--out_reads', required=True)
    parser.add_argument('--out_dupes', required=True)
    return parser.parse_args()


def get_duplicates(dupes_file):
    return [line.strip() for line in open(dupes_file)]


def get_reads(reads):
    return (line for line in SeqIO.parse(open(reads, 'r'), 'fastq'))


def make_out_file(fq):
    open(fq, 'w').close()


def write_out(outfile, read):
    with open(outfile , "a") as o:
        SeqIO.write(read, o, 'fastq')


def gen_deduped_reads(reads, duplicates, out_reads, out_dupes):
    for read in reads:
        if read.id in duplicates: write_out(out_dupes, read)
        else: write_out(out_reads, read)


def main():
    args = parse_args()
    duplicates = get_duplicates(args.duplicates)
    reads = get_reads(args.reads)
    make_out_file(args.out_reads)
    make_out_file(args.out_dupes)
    gen_deduped_reads(reads, duplicates, args.out_reads, args.out_dupes)


if __name__ == "__main__":
    main()