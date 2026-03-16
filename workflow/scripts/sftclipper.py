#!/usr/bin/env python3

from argparse import ArgumentParser
from itertools import groupby
import pysam

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--bam')
    parser.add_argument('--outfile')
    return parser.parse_args()


def parse_cigar(cigar):
    cigar_string = []
    cig_iter = groupby(cigar, lambda c: c.isdigit())
    for g, n in cig_iter:
        cigar_string.append((int("".join(n)), "".join(next(cig_iter)[1])))
    return cigar_string


# M	Match           no insertion or deletions, bases may not agree
# I	Insertion       additional base in query (not in reference)
# D	Deletion        query is missing base from reference
# =	Equal           no insertions or deletions, and bases agree
# X	Not Equal       no insertions or deletions, bases do not agree
# N	None            no query bases to align, an expected read gap (spliced read)
# S	Soft-Clipped    bases on end of read are not aligned but stored in SAM
# H	Hard-Clipped    bases on end of read are not aligned, not stored in SAM
# P	Padding	neither read nor reference has a base here


def clip_check(read, cutoff):
    cigar = parse_cigar(read.cigarstring)
    soft_clip = 0
    mapped = 0
    total = 0
    for c in cigar:
        total += c[0]
        if c[1] == 'S': soft_clip += c[0]
        if c[1] in 'MID=': mapped += c[0]
    return soft_clip/total < cutoff and mapped >= 100


def main():
    args = parse_args()
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    removed_name = args.outfile.split('.')[0]
    removed_bam = pysam.AlignmentFile(f'{removed_name}_REMOVED.bam', "wb", template=bamfile)
    with pysam.AlignmentFile(args.outfile, "wb", template=bamfile) as outf:
        [outf.write(read) if clip_check(read, args.cutoff) else removed_bam.write(read) for read in bamfile.fetch()]
    removed_bam.close()
    bamfile.close()



if __name__ == '__main__':
    main()