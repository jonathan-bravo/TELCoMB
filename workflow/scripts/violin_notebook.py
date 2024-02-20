#!/usr/bin/env python

from math import log10
import os.path
import sys
import configparser

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import SeqIO

import csv
import json
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--samples_list', nargs='+')
    parser.add_argument('--dedup_string')
    parser.add_argument('--config_file')
    parser.add_argument('--output_plot')
    return parser.parse_args()

# Functions
def read_megares_v2_ontology(config):
    # Create ontology dictionary from MEGARes ontology file
    megares_ontology = dict()
    hierarchy_dict = dict()
    with open(config['DATABASE']['MEGARES_ONTOLOGY'], 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv)
        for row in ontology_reader:
            # Skip column names
            if row[0] == "header":
                continue

            typ = row[1]
            cl = row[2]
            mech = row[3]
            group = row[4]

            # Set up hiearachy dict. This will be our tree structure
            if not typ in hierarchy_dict:
                hierarchy_dict[typ] = {}

            if not cl in hierarchy_dict[typ]:
                hierarchy_dict[typ][cl] = {}

            if not mech in hierarchy_dict[typ][cl]:
                hierarchy_dict[typ][cl][mech] = []

            if not group in hierarchy_dict[typ][cl][mech]:
                hierarchy_dict[typ][cl][mech].append(group)

            # FIll in our dict
            megares_ontology[row[0]] = {"class": cl,
                                        "mechanism": mech,
                                        "group": group,
                                        "type": typ
                                        }
    return megares_ontology, hierarchy_dict


def get_resistome(config, sample_name, dedup_string):
    return sample_name + dedup_string + config['EXTENSION']['RESISTOME_DIVERSITY']


def get_reads_length(config, sample_name):
    return sample_name + config['EXTENSION']['READS_LENGTH']


def main():
    args = parse_args()

    # Get info from snakemake
    samples_list = args.samples_list
    config_file = args.config_file
    config = configparser.ConfigParser()
    config.read(config_file)
    output_plot = args.output_plot

    # # Get info from snakemake
    # samples_list = snakemake.params[0]
    # config_file = snakemake.input[2]
    # config = configparser.ConfigParser()
    # config.read(config_file)
    # output_plot = snakemake.output[0]

    # Get megares genes size
    megares_gene_lengths = dict()
    for record in SeqIO.parse(config['DATABASE']['MEGARES'], "fasta"):
        megares_gene_lengths[record.description] = len(record.seq)

    # Get per sample size
    samples_size = dict()
    for sample in samples_list:
        samples_size[sample] = 0
        with open(get_reads_length(config, sample)) as json_file:
            reads_lengths = json.load(json_file)
            for read_name, length in reads_lengths.items():
                samples_size[sample] += length

    # Get per sample absolute abundance
    absolute_abundances = dict()

    for sample in samples_list:
        sample_absolute_abundance = dict()
        with open(get_resistome(config, sample, args.dedup_string)) as resistome_file:
            resistome_reader = csv.reader(resistome_file)
            for row in resistome_reader:
                # Remove statistics and headers
                if (row[0] == 'Statistics') or (row[0] == 'Resistome') or ('ARG_' in row[0]) or (row[0] == 'MEGARes Gene Header'):
                    continue

                gene_acc = row[0]
                gene_acc_hits = row[1]

                if 'RequiresSNPConfirmation' in gene_acc:
                    continue

                if gene_acc not in sample_absolute_abundance:
                    sample_absolute_abundance[gene_acc] = 0

                sample_absolute_abundance[gene_acc] += int(gene_acc_hits)

        absolute_abundances[sample] = sample_absolute_abundance

    # Convert absolute abundances to relative abundance
    relative_abundances = {}
    for sample, sample_absolute_abundance in absolute_abundances.items():
        sample_relative_abundance = {}
        for gene_acc, gene_acc_hits in sample_absolute_abundance.items():
            sample_relative_abundance[gene_acc] = log10((100 * float(gene_acc_hits)) / (float(megares_gene_lengths[gene_acc] * samples_size[sample])))

        relative_abundances[sample] = sample_relative_abundance
    
    # Plot
    df = pd.DataFrame.from_dict(relative_abundances)

    sns.set_style("whitegrid")
    sns.set_context("paper")
    ax = sns.violinplot(data=df, inner='box', palette='hls')
    ax.set(xlabel='Samples', ylabel='Log Relative Abundance')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=-30, ha="left")
    plt.gcf().subplots_adjust(bottom=0.30, right=0.85)
    ax.get_figure().savefig(output_plot)


if __name__ == '__main__':
    main()