#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--host_stats', required = True)
    parser.add_argument('--viral_stats', required = True)
    parser.add_argument('--outfile', required = True)
    return parser.parse_args()

def prep_host_stats(host_stats):
    df = pd.read_table(host_stats)
    df.drop(columns=['Unmapped'], inplace=True)
    df.columns = ['Sample', 'NumberOfInputReads', 'MappedToHost']
    return df

def prep_viral_stats(viral_stats):
    df = pd.read_table(viral_stats)
    df.drop(columns=['NumberOfInputReads', 'Unmapped'], inplace=True)
    df.columns=['Sample', 'MappedToViralDB']
    return df

def on_target_stats(host_df, viral_df, outfile):
    final_df = pd.merge(host_df, viral_df, on=['Sample'])
    final_df['Unmapped'] = final_df['NumberOfInputReads'] - (final_df['MappedToHost'] + final_df['MappedToViralDB'])
    final_df['HostReadsPercent'] = (final_df['MappedToHost']/final_df['NumberOfInputReads'])*100
    final_df['OnTargetPercent'] = (final_df['MappedToViralDB']/final_df['NumberOfInputReads'])*100
    final_df.to_csv(outfile, index=False, sep='\t')

def main():
    args = parse_args()
    host_df = prep_host_stats(args.host_stats)
    viral_df = prep_viral_stats(args.viral_stats)
    on_target_stats(host_df, viral_df, args.outfile)

if __name__ == '__main__':
    main()