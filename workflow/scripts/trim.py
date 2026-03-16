#!/usr/bin/env python3

import argparse
import gzip
from datetime import date
from Bio import SeqIO, bgzf
import re


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--crop', type=int, required=True)
    parser.add_argument('--logfile', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--barcodes', required=True)
    return parser.parse_args()


def get_reads(fq):
    return (seq for seq in SeqIO.parse(gzip.open(fq, 'rt'), 'fastq'))


def get_barcodes(fq):
    return [str(fasta.seq) for fasta in SeqIO.parse(open(fq), 'fasta')]


def make_log(crop): # change this to reflect changes to code
    with open(logfile, 'w') as o:
        o.write(f'{date.today()}\n\n')
        o.write(f'Base crop length is {crop} bases from head and tail of all reads\n')
        o.write('Some reads may have this value adjusted, this will be listed per read\n\n')


def make_fq():
    open(outfile, 'wb').close()


def write_seq_log(read, head_crop, tail_crop):
    with open(logfile, 'a') as o:
        o.write('sequence too short with crop lengths:\n')
        o.write(f'head crop len: {head_crop}, tail crop len: {-tail_crop}\n')
        o.write(f'>{read.id}\n')
        o.write(f'{str(read.seq)}\n\n')


def write_crop_log(read, head_crop, tail_crop):
    with open(logfile, 'a') as o:
        o.write(f'>{read.id}\n')
        o.write(f'head crop len: {head_crop}, tail crop len: {-tail_crop}\n')
        o.write(f'cropped head seq\n')
        o.write(f'{str(read.seq)[:head_crop]}\n')
        o.write(f'cropped tail seq\n')
        o.write(f'{str(read.seq)[tail_crop:]}\n\n')


def write_read(read, head_crop, tail_crop):
    with bgzf.BgzfWriter(outfile, 'ab') as o:
        SeqIO.write(read[head_crop:tail_crop], o, 'fastq')


def write_read_info():
    with open(logfile, 'a') as o:
        o.write('ReadID\tReadLength\tNumMatches\t(Head,Tail,Neither)\n')
        [o.write(f'{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}\n') for i in read_info]


def long_enough(seq, head_crop, tail_crop): # - tail_crop because it is negative
    return len(seq) >= (head_crop - tail_crop + 50)


def in_head(barcode):
    return 0 <= barcode[0] <= 50


def in_tail(barcode, seq_len):
    return seq_len - 50 <= barcode[0] <= seq_len


def find_matches(barcode, seq):
    # removed match.group() from tuple because it is the sequence of the barcode
    return [(m.start(), m.end()) for m in re.finditer(rf'{barcode}', seq)]


def barcode_check(seq):
    return [b for bar in barcodes for b in find_matches(bar, seq) if bar in seq]


def trim_read(read, head_crop, tail_crop):
    write_read(read, head_crop, tail_crop)
    write_crop_log(read, head_crop, tail_crop)


def point_read(read, head_crop, tail_crop):
    if long_enough(str(read.seq), head_crop, tail_crop):
        trim_read(read, head_crop, tail_crop)
    else:
        write_seq_log(read, head_crop, tail_crop)


def split_read(read, barcode):
    read_a = read[:barcode[0]]
    read_b = read[barcode[0]:]
    read_a.id = f'{read_a.id}_A'
    read_b.id = f'{read_b.id}_B'
    # read_a.description = read.description + f' new_read_id={read_a.id}_A'
    # read_b.description = read.description + f' new_read_id={read_b.id}_B'
    process_read(read_a)
    process_read(read_b)


def classify_match(match, seq_len):
    if in_head(match): return 'head'
    if in_tail(match, seq_len): return 'tail'
    return 'neither'


def bin_barcodes(matches, seq_len): # if this doesn't work I can use the below code

    bins = {'head': [], 'tail': [], 'neither': []}

    [bins[classify_match(m, seq_len)].append(m) for m in matches]

    head = bins['head']
    tail = bins['tail']
    neither = bins['neither']

    head.sort()
    tail.sort()
    neither.sort()

    return (head, tail, neither)


def process_read(read):
    head_crop = crop_len
    tail_crop = -crop_len
    
    matches = barcode_check(str(read.seq))

    head, tail, neither = bin_barcodes(matches, len(str(read.seq)))

    read_info.append((
        read.id,
        len(str(read.seq)),
        len(matches),
        (len(head), len(tail), len(neither))
    ))

    if len(neither) > 0:
        barcode = neither[0]
        split_read(read, barcode)
    else:
        if len(head) > 0: head_crop = head[-1][1] + 13
        if len(tail) > 0: tail_crop = tail[0][0] - 13 - len(str(read.seq))
        point_read(read, head_crop, tail_crop)


def main():
    # declaring globals
    global barcodes
    global crop_len
    global logfile
    global outfile
    global read_info

    # parsing input arguments
    args = parse_args()

    # getting all reads as generator object
    reads = get_reads(args.infile)

    # setting values for global variables
    barcodes = get_barcodes(args.barcodes)
    crop_len = args.crop
    logfile = args.logfile
    outfile = args.outfile
    read_info = []

    # make the inistial log and out file
    make_log(args.crop)
    make_fq()

    # begin processing
    [process_read(read) for read in reads]

    write_read_info()


if __name__ == '__main__':
    main()