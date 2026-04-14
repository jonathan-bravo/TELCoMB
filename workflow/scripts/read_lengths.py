#!/usr/bin/env python3

import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    return parser.parse_args()

def skip_lines(reads):
    next(reads)
    next(reads)

def process_fastq_read(read, reads):
    read_id = read.split(' ')[0][1:]
    read_len = len(next(reads))
    skip_lines(reads)
    return f'{read_id}\t{read_len}\n'

def get_read_lengths(reads):
    return [process_fastq_read(read, reads) for read in reads]

def write_out(outfile, read_lengths):
    with open(outfile, 'w') as o:
        o.write('ReadId\tReadLength\n')
        [o.write(length) for length in read_lengths]

def main():
    args = parse_args()
    reads = (line.strip() for line in gzip.open(args.infile, 'rt'))
    read_lengths = get_read_lengths(reads)
    write_out(args.outfile, read_lengths)

if __name__ == '__main__':
    main()