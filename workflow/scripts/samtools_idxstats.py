#!/usr/bin/env python3

import argparse
import re

def parse_cmdline_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', required=True)
    parser.add_argument('-o', required=True)
    return parser.parse_args()

def write_out(idxstats, output_file):
    with open(output_file, 'w') as o:
        o.write('Sample\tNumberOfInputReads\tMapped\tUnmapped\n')
        [o.write(stat) for stat in idxstats]

def mapping_stats(infile):
    # sample_name = infile.split('/')[-1].split('.', 1)[0]
    sample_name = re.split(r'_remove_host_samtools\.idxstats|_all_viruses_samtools\.idxstats', infile.split('/')[-1])[0]
    idxstats = ((int(row.strip().split('\t')[2]),int(row.strip().split('\t')[3])) for row in open(infile))
    mapped, unmapped = [sum(x) for x in zip(*idxstats)]
    number_of_reads = mapped + unmapped
    return f'{sample_name}\t{number_of_reads}\t{mapped}\t{unmapped}\n'

def main():
    args = parse_cmdline_params()
    idxstats = list(map(mapping_stats, args.i))
    write_out(idxstats, args.o)

if __name__ == "__main__":
    main()