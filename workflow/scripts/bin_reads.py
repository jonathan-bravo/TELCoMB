#!/usr/bin/env python

import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outdir', required=True)
    return parser.parse_args()

def gen_read_length_clusters(reads, outdir):
    for read in reads:
        read_id = read.split(' ')[0][1:]
        seq = next(reads)
        seq_len = len(seq.strip())
        next(reads) # desc
        next(reads) # qual
        outfile = f'{outdir}/{seq_len}.rl.bins.fasta.gz'
        with open(outfile, "a") as o:
            o.write('>' + read_id + '\n')
            o.write(seq)

def main():
    args = parse_args()
    reads = (line for line in open(args.infile, 'r'))
    gen_read_length_clusters(reads, args.outdir)

if __name__ == "__main__":
    main()