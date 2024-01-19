#!/usr/bin/env python

from Bio import SeqIO
import pysam
import argparse
import csv
import configparser
import numpy as np
import json
import os
from common import *


def long_reads_strategy_resistome(config):
    # Get megares lengths for coverage
    AMR_mapped_regions_per_read = dict()
    not_valid_AMR_mapped_region = dict()

    megares_gene_lengths = {}
    megares_reference_fasta_filename = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    sam_file = pysam.AlignmentFile(config['INPUT']['ARGS_SAM_FILE'], 'r')

    # Create ontology dictionary from MEGARes ontology file
    megares_ontology, _ = read_megares_v2_ontology(config)

    # Get reads lengths
    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_NO_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    # Iterate through every read. Accumulate number of reads while recording read length
    reads_aligned = set()
    genes_per_read = dict()
    for read in sam_file.fetch():
        if read.is_unmapped: continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'].upper() != 'TRUE' and read.is_secondary: continue

        # check coverage
        if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_AMR_THRESHOLD']):
            classname = megares_ontology[read.reference_name]["class"]
            mech = megares_ontology[read.reference_name]["mechanism"]
            group = megares_ontology[read.reference_name]["group"]

            # update gene dict
            if read.reference_name not in gene_dict:
                gene_dict[read.reference_name] = 1
            else:
                gene_dict[read.reference_name] += 1

            # update class dict
            if classname not in class_dict:
                class_dict[classname] = 1
            else:
                class_dict[classname] += 1

            # update mechanism dict
            if mech not in mech_dict:
                mech_dict[mech] = 1
            else:
                mech_dict[mech] += 1

            # update group dict
            if group not in group_dict:
                group_dict[group] = 1
            else:
                group_dict[group] += 1

            reads_aligned.add(read.query_name)
            if read.query_name not in genes_per_read:
                genes_per_read[read.query_name] = 0
            genes_per_read[read.query_name] += 1

            if read.query_name not in AMR_mapped_regions_per_read:
                AMR_mapped_regions_per_read[read.query_name] = list()
            AMR_mapped_regions_per_read[read.query_name].append([read.query_alignment_start, read.query_alignment_end, read.reference_name])

    # Mobilome
    logger = logging.getLogger()
    strategy_mob = config['MISC']['MOBILOME_STRATEGY']
    logger.info('Mobilome strategy: {}'.format(strategy_mob))
    if strategy_mob == 'SHORT':
        not_valid_AMR_mapped_region = short_reads_strategy_mobilome(config, AMR_mapped_regions_per_read)
    elif strategy_mob == 'LONG':
        not_valid_AMR_mapped_region = long_reads_strategy_mobilome(config, AMR_mapped_regions_per_read)
    else:
        logger.error("MOBILOME_STRATEGY must be in [LONG, SHORT], {} not supported".format(strategy_mob))

    for query_name in AMR_mapped_regions_per_read:
        for reference_info in AMR_mapped_regions_per_read[query_name]:
            if query_name in not_valid_AMR_mapped_region:
                if reference_info in not_valid_AMR_mapped_region[query_name]:
                    gene_dict[reference_info[2]] -= 1
                    if gene_dict[reference_info[2]] == 0:
                        gene_dict.pop(reference_info[2])
                    genes_per_read[query_name] -= 1
                    if genes_per_read[query_name] == 0:
                        reads_aligned.remove(read.query_name)
                    classname = megares_ontology[reference_info[2]]["class"]
                    class_dict[classname] -= 1
                    if class_dict[classname] == 0:
                        class_dict.pop(classname)
                    mech = megares_ontology[reference_info[2]]["mechanism"]
                    mech_dict[mech] -= 1
                    if mech_dict[mech] == 0:
                        mech_dict.pop(mech)
                    group = megares_ontology[reference_info[2]]["group"]
                    group_dict[group] -= 1
                    if group_dict[group] == 0:
                        group_dict.pop(group)
                    continue

    
    # Prepare rows of diversity csv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['ARG_' + stat_name, stat_value])

    csv_rows.append(['Resistome'])
    csv_rows.append(
        ["MEGARes Gene Header", "Num Reads", "Group", "Num Reads", "Mechanism", "Num Reads", "Class",
         "Num Reads"])
    start = len(csv_rows)
    for gene in sorted(gene_dict, key=lambda gene: gene_dict[gene], reverse=True):
        csv_rows.append([gene, gene_dict[gene]])

    i = start
    for group in sorted(group_dict, key=lambda group: group_dict[group], reverse=True):
        csv_rows[i].extend([group, group_dict[group]])
        i += 1

    i = start
    for mech in sorted(mech_dict, key=lambda mech: mech_dict[mech], reverse=True):
        csv_rows[i].extend([mech, mech_dict[mech]])
        i += 1

    i = start
    for class_name in sorted(class_dict, key=lambda class_name: class_dict[class_name], reverse=True):
        csv_rows[i].extend([class_name, class_dict[class_name]])
        i += 1

    # Write diversity tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_diversity.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    # Prepare rows of richness tsv
    csv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    csv_rows.append([len(gene_dict), len(class_dict), len(mech_dict), len(group_dict)])

    # Write richness tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_richness.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    return AMR_mapped_regions_per_read

def short_reads_strategy_resistome(config):
    # Get megares gene for coverage

    AMR_mapped_regions_per_read = dict()
    not_valid_AMR_mapped_region = dict()

    megares_genes = {}
    megares_reference_fasta_filename = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_genes[rec.name] = np.zeros(len(rec.seq))

    sam_file = pysam.AlignmentFile(config['INPUT']['ARGS_SAM_FILE'], 'r')

    # Create ontology dictionary from MEGARes ontology file
    megares_ontology, _ = read_megares_v2_ontology(config)

    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_NO_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    reads_aligned_per_gene = dict()
    # Iterate through every read. Accumulate number of reads while recording read length
    for read in sam_file.fetch():
        if read.is_unmapped: continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'].upper() != 'TRUE' and read.is_secondary: continue

        classname = megares_ontology[read.reference_name]["class"]
        mech = megares_ontology[read.reference_name]["mechanism"]
        group = megares_ontology[read.reference_name]["group"]

        # update gene dict
        if (not read.reference_name in gene_dict):
            gene_dict[read.reference_name] = 1
            reads_aligned_per_gene[read.reference_name] = dict()
            reads_aligned_per_gene[read.reference_name][read.query_name] = 1
        else:
            gene_dict[read.reference_name] += 1
            if read.query_name not in reads_aligned_per_gene[read.reference_name]:
                reads_aligned_per_gene[read.reference_name][read.query_name] = 0
            reads_aligned_per_gene[read.reference_name][read.query_name] += 1

        # update class dict
        if classname not in class_dict:
            class_dict[classname] = 1
        else:
            class_dict[classname] += 1

        # update mechanism dict
        if mech not in mech_dict:
            mech_dict[mech] = 1
        else:
            mech_dict[mech] += 1

        # update group dict
        if group not in group_dict:
            group_dict[group] = 1
        else:
            group_dict[group] += 1

        if read.query_name not in AMR_mapped_regions_per_read:
            AMR_mapped_regions_per_read[read.query_name] = list()
        AMR_mapped_regions_per_read[read.query_name].append([read.query_alignment_start, read.query_alignment_end, read.reference_name, read.reference_start, read.reference_end])

    # Mobilome
    logger = logging.getLogger()
    strategy_mob = config['MISC']['MOBILOME_STRATEGY']
    logger.info('Mobilome strategy: {}'.format(strategy_mob))
    if strategy_mob == 'SHORT':
        not_valid_AMR_mapped_region = short_reads_strategy_mobilome(config, AMR_mapped_regions_per_read)
    elif strategy_mob == 'LONG':
        not_valid_AMR_mapped_region = long_reads_strategy_mobilome(config, AMR_mapped_regions_per_read)
    else:
        logger.error("MOBILOME_STRATEGY must be in [LONG, SHORT], {} not supported".format(strategy_mob))

    for query_name in AMR_mapped_regions_per_read:
        for reference_info in AMR_mapped_regions_per_read[query_name]:
            if query_name in not_valid_AMR_mapped_region:
                if reference_info in not_valid_AMR_mapped_region[query_name]:
                    gene_dict[reference_info[2]] -= 1
                    reads_aligned_per_gene[reference_info[2]][query_name] -= 1
                    if reads_aligned_per_gene[reference_info[2]][query_name] == 0:
                        reads_aligned_per_gene[reference_info[2]].pop(query_name)
                    classname = megares_ontology[reference_info[2]]["class"]
                    class_dict[classname] -= 1
                    if class_dict[classname] == 0:
                        class_dict.pop(classname)
                    mech = megares_ontology[reference_info[2]]["mechanism"]
                    mech_dict[mech] -= 1
                    if mech_dict[mech] == 0:
                        mech_dict.pop(mech)
                    group = megares_ontology[reference_info[2]]["group"]
                    group_dict[group] -= 1
                    if group_dict[group] == 0:
                        group_dict.pop(group)
                    continue
            for i in range(reference_info[3], reference_info[4]):
                megares_genes[reference_info[2]][i] = 1

    # check coverage
    covered_genes = set()
    reads_aligned = set()
    for megares_gene, coverage_vector in megares_genes.items():
        if (float(sum(coverage_vector) / len(coverage_vector)) > float(config['MISC']['GLOBAL_AMR_THRESHOLD'])):
            covered_genes.add(megares_gene)
            reads_aligned.update(list(reads_aligned_per_gene[megares_gene].keys()))


    # Get only covered genes
    gene_covered_dict = dict()
    class_covered_dict = dict()
    mech_covered_dict = dict()
    group_covered_dict = dict()

    for gene in covered_genes:
        classname = megares_ontology[gene]["class"]
        mech = megares_ontology[gene]["mechanism"]
        group = megares_ontology[gene]["group"]

        gene_covered_dict[gene]         = gene_dict[gene]
        class_covered_dict[classname]   = class_dict[classname]
        mech_covered_dict[mech]         = mech_dict[mech]
        group_covered_dict[group]       = group_dict[group]


    # Prepare rows of diversity csv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['ARG_' + stat_name, stat_value])

    csv_rows.append(['Resistome'])
    csv_rows.append(
        ["MEGARes Gene Header", "Num Reads", "Group", "Num Reads", "Mechanism", "Num Reads", "Class",
         "Num Reads"])
    start = len(csv_rows)
    for gene in sorted(gene_covered_dict, key=lambda gene: gene_covered_dict[gene], reverse=True):
        csv_rows.append([gene, gene_covered_dict[gene]])

    i = start
    for group in sorted(group_covered_dict, key=lambda group: group_covered_dict[group], reverse=True):
        csv_rows[i].extend([group, group_covered_dict[group]])
        i += 1

    i = start
    for mech in sorted(mech_covered_dict, key=lambda mech: mech_covered_dict[mech], reverse=True):
        csv_rows[i].extend([mech, mech_covered_dict[mech]])
        i += 1

    i = start
    for class_name in sorted(class_covered_dict, key=lambda class_name: class_covered_dict[class_name], reverse=True):
        csv_rows[i].extend([class_name, class_covered_dict[class_name]])
        i += 1

    # Write diversity tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_diversity.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    # Prepare rows of richness tsv
    csv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    csv_rows.append([len(gene_covered_dict), len(class_covered_dict), len(mech_covered_dict), len(group_covered_dict)])

    # Write richness tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_richness.csv",'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    return AMR_mapped_regions_per_read


def long_reads_strategy_mobilome(config, AMR_mapped_regions_per_read):
    not_valid_AMR_mapped_region = dict()
    mges_gene_lengths = dict()
    mge_combined_reference_fasta_filename = config['DATABASE']['MGES']

    for rec in SeqIO.parse(mge_combined_reference_fasta_filename, "fasta"):
        mges_gene_lengths[rec.name] = len(rec.seq)

    # Get reads lengths
    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_NO_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    # Get list of overlapped MGEs
    overlapped_mges = list()
    with open(config['INPUT']['OVERLAP_LIST'], 'r') as overlap_list: 
        for mge in overlap_list:
            overlapped_mges.append(mge)

    reads_aligned_to_args = dict()
    if config['INPUT']['ARGS_SAM_FILE'] != '':
        args_gene_lengths = dict()
        args_reference_fasta_filename = config['DATABASE']['MEGARES']

        for rec in SeqIO.parse(args_reference_fasta_filename, "fasta"):
            args_gene_lengths[rec.name] = len(rec.seq)

        args_sam_file = pysam.AlignmentFile(config['INPUT']['ARGS_SAM_FILE'], 'r')

    mges_sam_file = pysam.AlignmentFile(config['INPUT']['MGES_SAM_FILE'], 'r')

    gene_dict = dict()
    reads_aligned = set()
    # Iterate through every aligned segment
    for read in mges_sam_file.fetch():
        if read.is_unmapped: continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'].upper() != 'TRUE' and read.is_secondary: continue

        if read.reference_name in overlapped_mges: continue

        overlapping_with_amr_gene = False
        current_length = 0 if read.query_name not in not_valid_AMR_mapped_region else len(not_valid_AMR_mapped_region[read.query_name])
        if read.query_name in AMR_mapped_regions_per_read:
            for region in AMR_mapped_regions_per_read[read.query_name]:
                if (region[1] < read.query_alignment_end) and (region[0] >= read.query_alignment_start) :
                    if read.query_name not in not_valid_AMR_mapped_region:
                        not_valid_AMR_mapped_region[read.query_name] = list()
                    not_valid_AMR_mapped_region[read.query_name].append(region[2])
                elif (region[1] == read.query_alignment_end) and (region[0] > read.query_alignment_start) :
                    if read.query_name not in not_valid_AMR_mapped_region:
                        not_valid_AMR_mapped_region[read.query_name] = list()
                    not_valid_AMR_mapped_region[read.query_name].append(region[2])
                elif (region[1] >= read.query_alignment_end) and (region[0] <= read.query_alignment_start) :
                    overlapping_with_amr_gene = True
                    break

        if overlapping_with_amr_gene and (read.query_name in not_valid_AMR_mapped_region):
            while len(not_valid_AMR_mapped_region[read.query_name]) > current_length:
                not_valid_AMR_mapped_region[read.query_name].pop(current_length)
            continue

        # check coverage
        if (float(read.reference_length) / mges_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_MGE_THRESHOLD']):
            if not read.reference_name in gene_dict:
                gene_dict[read.reference_name] = 1
            else:
                gene_dict[read.reference_name] += 1
            reads_aligned.add(read.query_name)

    # Prepare rows of tsv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['MGES_' + stat_name, stat_value])

    csv_rows.append(['Mobilome'])

    gene_riches = [(header, gene_dict[header]) for header in gene_dict]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(gene_riches)])

    # Column headers
    csv_rows.append(["MGE Header", "Num Reads"])

    # Output individual MGEs and their conts
    for gene_count_tuple in sorted(gene_riches, key=lambda gene_count_tuple: gene_count_tuple[1], reverse=True):
        csv_rows.append([gene_count_tuple[0], gene_count_tuple[1]])

    # Write csv
    filename_prefix = config['OUTPUT']['OUTPUT_PREFIX']
    with open(filename_prefix + '_' + config['MISC']['MOBILOME_STRATEGY'] + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)
    return not_valid_AMR_mapped_region
    

def short_reads_strategy_mobilome(config, AMR_mapped_regions_per_read):
    not_valid_AMR_mapped_region = dict()
    mge_genes = dict()
    mge_combined_reference_fasta_filename = config['DATABASE']['MGES']
    for rec in SeqIO.parse(mge_combined_reference_fasta_filename, "fasta"):
        mge_genes[rec.name] = np.zeros(len(rec.seq))

    mges_sam_file = pysam.AlignmentFile(config['INPUT']['MGES_SAM_FILE'], 'r')

    # Get reads lengths
    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_NO_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    # Get list of overlapped MGEs
    overlapped_mges = list()
    with open(config['INPUT']['OVERLAP_LIST'], 'r') as overlap_list: 
        for mge in overlap_list:
            overlapped_mges.append(mge[:-1])

    gene_hits = dict()
    reads_aligned_per_gene = dict()
    # Iterate through every read. Accumulate number of reads aligned and number of alignments per aclame mge
    for read in mges_sam_file.fetch():
        if read.is_unmapped: continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'].upper() != 'TRUE' and read.is_secondary: continue

        if read.reference_name in overlapped_mges: continue

        overlapping_with_amr_gene = False
        current_length = 0 if read.query_name not in not_valid_AMR_mapped_region else len(not_valid_AMR_mapped_region[read.query_name])
        if read.query_name in AMR_mapped_regions_per_read:
            for region in AMR_mapped_regions_per_read[read.query_name]:
                if (region[1] < read.query_alignment_end) and (region[0] >= read.query_alignment_start) :
                    if read.query_name not in not_valid_AMR_mapped_region:
                        not_valid_AMR_mapped_region[read.query_name] = list()
                    not_valid_AMR_mapped_region[read.query_name].append(region[2])
                elif (region[1] == read.query_alignment_end) and (region[0] > read.query_alignment_start) :
                    if read.query_name not in not_valid_AMR_mapped_region:
                        not_valid_AMR_mapped_region[read.query_name] = list()
                    not_valid_AMR_mapped_region[read.query_name].append(region[2])
                elif (region[1] >= read.query_alignment_end) and (region[0] <= read.query_alignment_start) :
                    overlapping_with_amr_gene = True
                    break

        if overlapping_with_amr_gene and (read.query_name in not_valid_AMR_mapped_region):
            while len(not_valid_AMR_mapped_region[read.query_name]) > current_length:
                not_valid_AMR_mapped_region[read.query_name].pop(current_length)
            continue

        for i in range(read.reference_start, read.reference_end):
            mge_genes[read.reference_name][i] = 1

        if read.reference_name not in gene_hits:
            gene_hits[read.reference_name] = 1
            reads_aligned_per_gene[read.reference_name] = set()
            reads_aligned_per_gene[read.reference_name].add(read.query_name)
        else:
            gene_hits[read.reference_name] += 1
            reads_aligned_per_gene[read.reference_name].add(read.query_name)

    # check coverage
    covered_genes = set()
    reads_aligned = set()
    for mge_gene, coverage_vector in mge_genes.items():
        if float(sum(coverage_vector) / len(coverage_vector)) > float(config['MISC']['GLOBAL_MGE_THRESHOLD']):
            covered_genes.add(mge_gene)
            reads_aligned.update(reads_aligned_per_gene[mge_gene])

    # Prepare rows of tsv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['MGES_' + stat_name, stat_value])

    csv_rows.append(['Mobilome'])
    covered_gene_richness = [(header, gene_hits[header]) for header in covered_genes]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(covered_gene_richness)])

    # Column headers
    csv_rows.append(["MGE Header", "Num Reads"])

    # Output individual MGEs and their counts
    for gene_count_tuple in sorted(covered_gene_richness, key=lambda gene_count_tuple: gene_count_tuple[1], reverse=True):
        csv_rows.append([gene_count_tuple[0], gene_count_tuple[1]])

    # Write csv
    filename_prefix = config['OUTPUT']['OUTPUT_PREFIX']
    with open(filename_prefix + '_' + config['MISC']['MOBILOME_STRATEGY'] + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)
    return not_valid_AMR_mapped_region


def main():
    parser = argparse.ArgumentParser(description='Compute resistome')
    parser.add_argument('-r', help='Reads file', dest='reads_file', required=True)
    parser.add_argument('-a', help='ARGS Alignment file', dest='args_sam', required=True)
    parser.add_argument('-m', help='MGES Alignment file', dest='mges_sam', required=True)
    parser.add_argument('-o', help='Output Prefix', dest='out_prefix', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    parser.add_argument('-s', help='Overlapped MGEs list', dest='overlap', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    logger = init_logger()

    config['INPUT'] = dict()
    config['INPUT']['ARGS_SAM_FILE'] = args.args_sam
    config['INPUT']['MGES_SAM_FILE'] = args.mges_sam
    config['INPUT']['OVERLAP_LIST'] = args.overlap
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(args.reads_file)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(args.reads_file))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])

    config['OUTPUT'] = dict()
    config['OUTPUT']['OUTPUT_PREFIX'] = args.out_prefix
    config['OUTPUT']['OUT_DIR'] = os.path.dirname(os.path.abspath(config['OUTPUT']['OUTPUT_PREFIX']))

    # Resistome
    strategy_res = config['MISC']['RESISTOME_STRATEGY']
    logger.info('Resistome strategy: {}'.format(strategy_res))
    if strategy_res == 'SHORT':
        AMR_mapped_regions_per_read = short_reads_strategy_resistome(config)
    elif strategy_res == 'LONG':
        AMR_mapped_regions_per_read = long_reads_strategy_resistome(config)
    else:
        logger.error("RESISTOME_STRATEGY must be in [LONG, SHORT], {} not supported".format(strategy_res))




if __name__ == "__main__":
    main()

