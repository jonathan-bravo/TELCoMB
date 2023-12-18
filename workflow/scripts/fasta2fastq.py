#!/usr/bin/env python

from argparse import ArgumentParser
from Bio import SeqIO


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--fasta')
    parser.add_argument('--fastq')
    return parser.parse_args()


def read_fasta(fasta):
    return ((record.id, str(record.seq)) for record in SeqIO.parse(open(fasta, 'r'), "fasta"))

def write_fastq(fasta, fastq):
    with open(fastq, 'w') as f:
        for record in fasta:
            f.write(f'@{record[0]}\n')
            f.write(f'{record[1]}\n')
            f.write('+\n')
            f.write('~' * len(record[1]) + '\n')


def main():
    args = parse_args()
    fasta = read_fasta(args.fasta)
    write_fastq(fasta, args.fastq)


if __name__ == '__main__':
    main()