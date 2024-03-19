#!/usr/bin/env python

import argparse
from math import ceil
from os import listdir, system

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', required = True)
    parser.add_argument('--indir', required = True)
    parser.add_argument('--outdir', required = True)
    return parser.parse_args()

def get_bins(indir):
    bins = [int(f.split('.')[0]) for f in listdir(indir) if f != '.DS_Store']
    bins.sort()
    return bins

def get_clusters(bins, threshold):
    val = 1.99 - threshold
    clusters = []
    while bins:
        cluster_len = ceil(bins[0]*val)
        cluster = [x for x in bins if x <= cluster_len]
        clusters.append(cluster)
        bins = [x for x in bins if x not in cluster]
    return clusters

def cat_files(clusters, indir, outdir):
    for cluster in clusters:
        start = cluster[0]
        end = cluster[-1]
        files = ' '.join([f'{indir}/{x}.rl.bins.fasta.gz' for x in cluster])
        system(f'cat {files} > {outdir}/{start}.to.{end}.rl.clusters.fasta.gz')

def main():
    args = parse_args()
    bins = get_bins(args.indir)
    clusters = get_clusters(bins, args.threshold)
    cat_files(clusters, args.indir, args.outdir)

if __name__ == '__main__':
    main()