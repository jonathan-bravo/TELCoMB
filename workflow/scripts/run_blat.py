#!/usr/bin/env python3

from argparse import ArgumentParser
from os import listdir, getcwd, system
from concurrent.futures import ProcessPoolExecutor as ppe

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--threads', type = int, required = True)
    parser.add_argument('--outdir', required = True)
    parser.add_argument('--read_clusters', required = True)
    return parser.parse_args()

def get_read_clusters(rc):
    return [f'{getcwd()}/{rc}/{f}' for f in listdir(rc)]

def get_out_list(outdir, l):
    return [outdir] * l

def run_blat(cluster, outdir):
    cluster_id = cluster.split('/')[-1].split('_rl_clusters.fasta')[0]
    system(f'pblat -threads=10 -fastMap {cluster} {cluster} {outdir}{cluster_id}.psl')

def thread(threads, rc, out_list):
    with ppe(threads) as p:
        p.map(run_blat, rc, out_list)

def main():
    args = parse_args()
    rc = get_read_clusters(args.read_clusters)
    out_list = get_out_list(args.outdir, len(rc))
    thread(args.threads, rc, out_list)

if __name__ == '__main__':
     main()
