#!/usr/bin/env python

# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from common import *
import os


def main():
    parser = argparse.ArgumentParser(description='Plot read lengths')
    parser.add_argument('-s', help='Reject outliers. Standard Deviations', default=0, dest='std_dev', type=int)
    parser.add_argument('-b', help='Number of bins', default=100, dest='bins', type=int)
    parser.add_argument('-o', help="Output file", default='', dest='out_file', type=str)
    parser.add_argument('--title', help='Plot title', default='', dest='title', type=str)
    parser.add_argument('--inset', help="Add inset. Left endpoint of the interval", dest='inset', default=0, type=int)
    parser.add_argument('--inset-position', help='Inset position', default='0.3,0.5,0.6,0.4', dest='inset_pos', type=str)
    parser.add_argument('--x-scale', help='Set X axis scale, pyplot argument [log,linear]', dest='x_scale', default='linear')
    parser.add_argument('--y-scale', help='Set Y axis scale, pyplot argument [log,linear]', dest='y_scale', default='linear')
    parser.add_argument('-i', '--input', nargs='+', help='Input Fastq files.', dest='input_files', required=True)
    args = parser.parse_args()

    # Sanity check on input files
    for file_path in args.input_files:
        if not os.path.exists(file_path):
            print('The file {} does not exist!'.format(file_path))
            exit(FileNotFoundError)

    read_lengths = list()
    inset_read_lengths = list()
    for file_path in args.input_files:
        if is_gz_file(file_path):
            handle = gzip.open(file_path, 'rt')
        else:
            handle = open(file_path, 'rt')

        for record in SeqIO.parse(handle, "fastq"):
            read_length = len(record.seq)
            read_lengths.append(read_length)
            if args.inset != 0 and read_length > args.inset:
                inset_read_lengths.append(read_length)

    if args.std_dev != 0:
        read_lengths = reject_outliers(read_lengths, args.std_dev)
        max_read_length = max(read_lengths)
        inset_read_lengths = [i_read_length for i_read_length in inset_read_lengths if i_read_length <= max_read_length]

    # Plot histogram
    fig, ax = plt.subplots()
    ax.set_title(args.title)
    ax.set_xlabel("Read length (nt)")
    ax.set_ylabel("Frequency (Number of Reads)")
    ax.hist(read_lengths, bins=args.bins)

    # Set Axis scale
    ax.set_xscale(args.x_scale)
    ax.set_yscale(args.y_scale)

    if args.inset != 0 and len(inset_read_lengths) > 0:
        # Parse position arguments
        # Example default '0.3,0.5,0.6,0.4' translates to 
        # 60% width, 40% height, positioned at (0.3, 0.5)
        pos = [float(s) for s in args.inset_pos.split(',')]
        
        # Create the inset axes directly on 'ax'
        inset_ax = inset_axes(ax, 
                              width=f"{pos[2]*100}%", 
                              height=f"{pos[3]*100}%", 
                              loc='lower left',
                              bbox_to_anchor=(pos[0], pos[1], 1, 1),
                              bbox_transform=ax.transAxes)

        # Plot the histogram in the inset
        inset_ax.hist(inset_read_lengths, bins=args.bins)

    # Save
    if args.out_file != '':
        plt.savefig(args.out_file)
    else:
        plt.savefig(args.input_files[0] + '.hist.pdf')


if __name__ == '__main__':
    main()
