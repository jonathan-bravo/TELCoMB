#!/usr/bin/env python3

import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outdir', required=True)
    return parser.parse_args()

def gen_read_length_clusters(reads, outdir):
    for read in reads:
        read_id = read.split(b' ')[0][1:]
        seq = next(reads)
        seq_len = len(seq.strip())
        next(reads) # desc
        next(reads) # qual
        outfile = f'{outdir}/{seq_len}_rl_bins.fasta'
        with open(outfile, "ab") as o:
            o.write(b'>' + read_id + b'\n')
            o.write(seq)

def main():
    args = parse_args()
    reads = (line for line in gzip.open(args.infile, 'rb'))
    gen_read_length_clusters(reads, args.outdir)

if __name__ == "__main__":
    main()
