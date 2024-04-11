#!/usr/bin/env python

import seaborn as sns
import csv
import json
import configparser
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--config_file')
    parser.add_argument('--read_lengths')
    parser.add_argument('--colocalizations')
    parser.add_argument('--output_plot')
    return parser.parse_args()


# def read_megares_v2_ontology(config):
#     # Create ontology dictionary from MEGARes ontology file
#     megares_ontology = dict()
#     hierarchy_dict = dict()
#     with open(config['DATABASE']['MEGARES_ONTOLOGY'], 'r') as ontology_tsv:
#         ontology_reader = csv.reader(ontology_tsv)
#         for row in ontology_reader:
#             # Skip column names
#             if row[0] == "header":
#                 continue

#             typ = row[1]
#             cl = row[2]
#             mech = row[3]
#             group = row[4]

#             # Set up hiearachy dict. This will be our tree structure
#             if not typ in hierarchy_dict:
#                 hierarchy_dict[typ] = {}

#             if not cl in hierarchy_dict[typ]:
#                 hierarchy_dict[typ][cl] = {}

#             if not mech in hierarchy_dict[typ][cl]:
#                 hierarchy_dict[typ][cl][mech] = []

#             if not group in hierarchy_dict[typ][cl][mech]:
#                 hierarchy_dict[typ][cl][mech].append(group)

#             # FIll in our dict
#             megares_ontology[row[0]] = {"class": cl,
#                                         "mechanism": mech,
#                                         "group": group,
#                                         "type": typ
#                                         }
#     return megares_ontology, hierarchy_dict
    

def add_genes_to_plot(gene_list, color, label, axis):
    for gene, position in gene_list:
        length = int(position.split(':')[1]) - int(position.split(':')[0])
        start = int(position.split(':')[0])
        if label == 'ARG':
            gene_label = gene.split('|')[4]
        elif label == 'MGE':
            gene_label = gene
        element_r = Rectangle((start, curr_y), length, height, facecolor=color, alpha=0.8, edgecolor='black',label=label)
        axis.add_patch(element_r)

        rx, ry = element_r.get_xy()
        cx = rx + element_r.get_width() / 2.0
        cy = ry + element_r.get_height() / 2.0
        axis.annotate(gene_label, (cx, cy - (height / 2 + 1)), color='black', fontsize=8, ha='center',va='center')


def main():
    args = parse_args()

    # Get info from snakemake
    config_file = args.config_file
    config = configparser.ConfigParser()
    config.read(config_file)
    reads_lengths = args.read_lengths
    colocalizations = args.colocalizations
    output_plot = args.output_plot
    # max_length = config["MISC"]["MAX_BP_COLOCALIZATIONS_PLOT"]

    global curr_y
    global height

    # Get read lengths
    with open(reads_lengths, 'r') as json_file:
        reads_lengths = json.load(json_file)

    ARG_COLOR  = 'red'
    MGE_COLOR  = 'green'
    # READ_COLOR = 'blue'

    # megares_ontology, _ = read_megares_v2_ontology(config)

    colocalizations_list = [row for row in csv.reader(open(colocalizations, 'rt'))][1:]

    used_reads_lengths = list()
    padding = 3
    height = 2
    curr_y = 1 + padding
    plt.figure(figsize=(16, 3*(len(colocalizations_list)/7)))
    currentAxis = plt.gca()
    for row in colocalizations_list:
        read_name = row[0]
        read_length = reads_lengths[read_name]
        used_reads_lengths.append(read_length)
        # read_r = Rectangle((0, curr_y), read_length, height, facecolor='blue', alpha=0.8, edgecolor='black', label='Read')
        read_r = Rectangle((0, curr_y+0.75), read_length, 0.5, facecolor='black', alpha=0.8, edgecolor='black', label='Read')
        currentAxis.add_patch(read_r)

        # Get ARGs, MGEs, KEGGs
        args = zip(row[1].split(';'), row[2].split(';'))
        mges = zip(row[3].split(';'), row[4].split(';'))

        # Add ARGs to plot
        add_genes_to_plot(args, ARG_COLOR, 'ARG', currentAxis)

        # Add MGEs to plot
        add_genes_to_plot(mges, MGE_COLOR, 'MGE', currentAxis)

        curr_y += height + padding


    sns.despine(top=True, right=True, left=True, bottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    plt.xlim([0, max(used_reads_lengths)])
    plt.ylim([0, curr_y])
    # plt.grid(axis='x', color='0.95')
    plt.xlabel('Read')

    # remove duplicates from legends
    handles, labels = currentAxis.get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)

    plt.legend(newHandles, newLabels)
    plt.savefig(output_plot)


if __name__ == '__main__':
    main()