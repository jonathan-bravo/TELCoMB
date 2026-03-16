#!/usr/bin/env python3

from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--run_id', required=True)
    parser.add_argument('--infiles', nargs='+', required=True)
    return parser.parse_args()


def rewrite_file(infile, run_id, sample_id):
    original = [line for line in open(infile,'r')]
    with open(infile, 'w') as final:
        final.write(f'run_id\tsample_id\t{original[0]}')
        for l in original[1:]:
            final.write(f'{run_id}\t{sample_id}\t{l}')


def main():
    args = parse_args()
    for f in args.infiles:
        rewrite_file(f, args.run_id, args.sample_id)


if __name__ == '__main__':
    main()