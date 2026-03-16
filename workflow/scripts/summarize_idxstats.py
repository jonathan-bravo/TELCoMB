#!/usr/bin/env python3

import argparse

def parse_cmdline_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--strains', required=True)
    parser.add_argument('--outfile', required=True)
    return parser.parse_args()

def parse_strains(strain_db):
    strains = (row.split('\t') for row in open(strain_db))
    return dict([(row[0], row[1]) for row in strains])

def get_strain_name(strains, strain):
    try: name = strains[strain]
    except KeyError: name = '*'
    return name.strip()

def write_out(idxstats, strains, output_file):
    with open(output_file, 'w') as o:
        # o.write('Strain\tName\tNumberOfInputReads\tMapped\tUnmapped\tPercentMapped\n')
        # [o.write(f'{strain}\t{get_strain_name(strains, strain)}\t{idxstats[strain][0]}\t{idxstats[strain][1]}\t{idxstats[strain][2]}\t{idxstats[strain][3]}\n') for strain in idxstats]
        o.write('Strain\tName\tNumberOfMappedReads\n')
        [o.write(f'{strain}\t{get_strain_name(strains, strain)}\t{idxstats[strain][0]}\n') for strain in idxstats]


def mapping_stats(infile):
    idxstats = (row.strip().split('\t') for row in open(infile))
    stats = {}
    for line in idxstats:
        number_of_reads = int(line[2]) + int(line[3])
        l1 = [number_of_reads, int(line[2]), int(line[3])]
        try:
            l2 = stats[line[0]]
            stats[line[0]] = [sum(x) for x in zip(*[l1,l2])]
        except KeyError:
            stats[line[0]] = l1
    return stats

def main():
    args = parse_cmdline_params()
    strains = parse_strains(args.strains)
    idxstats = mapping_stats(args.infile) # dictionary
    [idxstats[strain].append((idxstats[strain][1]/idxstats[strain][0])*100) if idxstats[strain][0] != 0 else idxstats[strain].append(0) for strain in idxstats] # Add name here?
    write_out(idxstats, strains, args.outfile)

if __name__ == "__main__":
    main()