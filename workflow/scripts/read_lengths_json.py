#!/usr/bin/env python

# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

import json
from Bio import SeqIO
from common import *

def main():
    parser = argparse.ArgumentParser(description='Get reads lengths from reads file')
    parser.add_argument('input', help='Input Fastq file')
    args = parser.parse_args()

    file_path = args.input
    if is_gz_file(file_path):
        handle = gzip.open(file_path, 'rt')
    else:
        handle = open(file_path, 'rt')

    read_lengths = dict()
    for record in SeqIO.parse(handle, "fastq"):
        read_lengths[record.name] = len(record.seq)

    print(json.dumps(read_lengths))


if __name__ == '__main__':
    main()
