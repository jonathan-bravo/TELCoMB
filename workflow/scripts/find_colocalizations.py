#!/usr/bin/env python

from Bio import SeqIO
import pysam
import argparse
from common import *


# Check if valid alignements, expecting intervals to be sorted
def not_overlapping(intervals_list):
    sorted_intervals = sorted(intervals_list)

    curr_right_lim = 0
    for interval in sorted_intervals:
        if interval[0] < curr_right_lim:
            return False
        curr_right_lim = interval[1]
    return True


# def get_colocalizations(config, reads_file_path, to_megares_path, to_mges_path, to_kegg_path):
def get_colocalizations(config, reads_file_path, to_megares_path, to_mges_path):
    logger = logging.getLogger()

    # Get read lengths
    logger.info("Reading FASTQ")
    reads_length = dict()
    if is_gz_file(reads_file_path):
        reads_file_handle = gzip.open(reads_file_path, 'rt')
    else:
        reads_file_handle = open(reads_file_path, 'rt')

    for record in SeqIO.parse(reads_file_handle, "fastq"):
        reads_length[record.name] = len(record.seq)

    reads_file_handle.close()

    overlapped_mges = list()
    with open(config['INPUT']['OVERLAP_LIST'], 'r') as overlap_list: 
        for mge in overlap_list:
            overlapped_mges.append(mge[:-1])

    # Get AMR genes lengths for coverage
    logger.info("Reading ARGS DB")
    megares_gene_lengths = {}
    megares_reference_fasta_filename = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    # Get MGEs lengths for coverage
    logger.info("Reading MGES DB")
    mge_gene_lengths = dict()
    mges_reference_fasta_filename = config['DATABASE']['MGES']
    for rec in SeqIO.parse(mges_reference_fasta_filename, "fasta"):
        mge_gene_lengths[rec.name] = len(rec.seq)

    # Get KEGG lengths for coverage
    # kegg_gene_lengths = dict()
    # kegg_reference_fasta_filename = config['DATABASE']['KEGG']
    # for rec in SeqIO.parse(kegg_reference_fasta_filename, "fasta"):
    #     kegg_gene_lengths[rec.name] = len(rec.seq)

    # Open aligned to megares sam
    logger.info("Reading ARGS alignment file")
    amr_positions = dict()
    read_to_amr = dict()
    amr_to_generated_bases = dict()
    amr_length = dict()
    with pysam.AlignmentFile(to_megares_path, "r") as to_megares_samfile:
        for read in to_megares_samfile:
            if read.is_unmapped: continue

            if config['MISC']['USE_SECONDARY_ALIGNMENTS'].upper() != 'TRUE' and read.is_secondary: continue
            # if read.is_secondary:
            #     continue

            # Check coverage
            if (read.reference_length / (megares_gene_lengths[read.reference_name])) > float(config['MISC']['GLOBAL_AMR_THRESHOLD_COLOCALIZATIONS']):
                if read.query_name not in read_to_amr:
                    read_to_amr[read.query_name] = list()
                    amr_positions[read.query_name] = list()
                amr_positions[read.query_name].append([read.query_alignment_start, read.query_alignment_end])
                read_to_amr[read.query_name].append(read.reference_name)

                amr_length[read.reference_name] = megares_gene_lengths[read.reference_name]

                if read.reference_name not in amr_to_generated_bases:
                    amr_to_generated_bases[read.reference_name] = read.reference_length
                else:
                    amr_to_generated_bases[read.reference_name] += read.reference_length

    # Open aligned to MGEs
    logger.info("Reading MGEs alignment files")
    mge_positions = dict()
    read_to_mges = dict()
    with pysam.AlignmentFile(to_mges_path, "r") as mge_alignment_file:
        for read in mge_alignment_file:
            if read.is_unmapped: continue

            if config['MISC']['USE_SECONDARY_ALIGNMENTS'].upper() != 'TRUE' and read.is_secondary: continue
            # if read.is_secondary:
            #     continue

            if read.reference_name in overlapped_mges: continue

            # Check coverage
            gene_length = mge_gene_lengths[read.reference_name]

            if (read.reference_length / gene_length) > float(config['MISC']['GLOBAL_MGE_THRESHOLD_COLOCALIZATIONS']):
                if read.query_name not in read_to_mges:
                    read_to_mges[read.query_name] = list()
                    mge_positions[read.query_name] = list()
                mge_positions[read.query_name].append(
                    [read.query_alignment_start, read.query_alignment_end])
                read_to_mges[read.query_name].append(read.reference_name)

    genes_lists = dict()
    genes_list_csv = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION'][
        'GENES_LIST']
    logger.info("Writing per read genes list to {}".format(genes_list_csv))
    with open(genes_list_csv, 'w') as genes_list_handle:
        writer = csv.writer(genes_list_handle)
        # header = ['Read Name',
        #           'AMR Genes', 'AMR Genes Pos',
        #           'MGE Genes', 'MGE Genes Pos',
        #           'KEGG genes', 'KEGG Genes Pos'
        #           ]
        header = ['Read Name',
                  'AMR Genes', 'AMR Genes Pos',
                  'MGE Genes', 'MGE Genes Pos'
                  ]
        writer.writerow(header)

        for read, amr_list in read_to_amr.items():
            # Amr genes list
            amr_genes_list = list()
            for idx, amr_name in enumerate(amr_list):
                amr_genes_list.append([amr_name, amr_positions[read][idx], 'amr'])

            # # Kegg genes list
            # kegg_genes_list = list()
            # if read in read_to_kegg:
            #     for idx, kegg_name in enumerate(read_to_kegg[read]):
            #         kegg_genes_list.append([kegg_name, kegg_positions[read][idx], 'kegg'])

            # MGE genes list
            mge_genes_list = list()
            if read in read_to_mges:
                for idx, mge_name in enumerate(read_to_mges[read]):
                    mge_genes_list.append([mge_name, mge_positions[read][idx], 'mge'])

            genes_lists[read] = list()
            genes_lists[read].extend(amr_genes_list)
            genes_lists[read].extend(mge_genes_list)
            # genes_lists[read].extend(kegg_genes_list)

            amr_genes_concatenated_names = list()
            amr_genes_concatenated_pos = list()
            for name, pos, _ in amr_genes_list:
                amr_genes_concatenated_names.append(name)
                amr_genes_concatenated_pos.append('({},{})'.format(pos[0], pos[1]))

            mge_genes_concatenated_names = list()
            mge_genes_concatenated_pos = list()
            for name, pos, _ in mge_genes_list:
                mge_genes_concatenated_names.append(name)
                mge_genes_concatenated_pos.append('({},{})'.format(pos[0], pos[1]))

            # kegg_genes_concatenated_names = list()
            # kegg_genes_concatenated_pos = list()
            # for name, pos, _ in kegg_genes_list:
            #     kegg_genes_concatenated_names.append(name)
            #     kegg_genes_concatenated_pos.append('({},{})'.format(pos[0], pos[1]))

            # row = [read,
            #        ';'.join(amr_genes_concatenated_names), ';'.join(amr_genes_concatenated_pos),
            #        ';'.join(mge_genes_concatenated_names), ';'.join(mge_genes_concatenated_pos),
            #        ';'.join(kegg_genes_concatenated_names), ';'.join(kegg_genes_concatenated_pos)]
            row = [read,
                   ';'.join(amr_genes_concatenated_names), ';'.join(amr_genes_concatenated_pos),
                   ';'.join(mge_genes_concatenated_names), ';'.join(mge_genes_concatenated_pos)]
            writer.writerow(row)

    # Candidate colocalizations
    candidate_colocalizations = {k: v for k, v in genes_lists.items() if len(v) >= 2}

    colocalizations = dict()
    for read, candidate_colocalization_list in candidate_colocalizations.items():
        if read not in reads_length:
            logger.error("{} not in read lengths.".format(read))
            continue

        has_ARG = False
        has_MGE = False
        colocalization = list()
        sorted_candidate_coloc_list = sorted(candidate_colocalization_list, key=lambda x: (x[1], x[0]))
        #sorted_candidate_coloc_list = sorted(candidate_colocalization_list, key=lambda x: x[0])

        # First take all the non-overlapping ARGS
        for idx, aligned_gene in enumerate(sorted_candidate_coloc_list):
            # if aligned_gene[2] == 'kegg':
            #     continue
            if len(colocalization) == 0 or aligned_gene[1][0] > colocalization[-1][1][1]:
                colocalization.append(aligned_gene)
                if aligned_gene[2] == 'amr':
                    has_ARG = True
                if aligned_gene[2] == 'mge':
                    has_MGE = True
            elif (aligned_gene[1][0] > colocalization[-1][1][0]) and (aligned_gene[1][1] > colocalization[-1][1][1]):
                aligned_gene[1][0] = colocalization[-1][1][1]+1
                colocalization.append(aligned_gene)
                if aligned_gene[2] == 'amr':
                    has_ARG = True
                if aligned_gene[2] == 'mge':
                    has_MGE = True

        # Finally fit all non overlapping genes from Kegg
        # for idx, aligned_gene in enumerate(sorted_candidate_coloc_list):
        #     if aligned_gene[2] == 'kegg':
        #         ti_coloc = [c[1] for c in colocalization]
        #         ti_coloc.append(aligned_gene[1])
        #         if not_overlapping(ti_coloc):
        #             colocalization.append(aligned_gene)

        if has_ARG and has_MGE:
            colocalizations[read] = colocalization

    arg_set = set()
    for read, args_list in read_to_amr.items():
        for arg in args_list:
            arg_set.add(arg)

    mge_set = set()
    for read, mges_list in read_to_mges.items():
        for mge in mges_list:
            mge_set.add(mge)

    logger.info("MGES: {}\tARGS: {}".format(len(mge_set), len(arg_set)))
    return arg_set, mge_set, colocalizations


############################################################
# Compute values

def main():
    parser = argparse.ArgumentParser(description='Colocalizations Finder.')
    # parser.add_argument('-k', help='Kegg alignment file', dest='kegg_sam', required=True)
    parser.add_argument('--arg', help='Megares alignment file', dest='megares_sam', required=True)
    parser.add_argument('--mge', help='MGEs alignment file', dest='mges_sam', required=True)
    parser.add_argument('-r', help='Reads file', dest='reads_file', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    parser.add_argument('-o', help='Output directory', dest='output_dir_path', required=True)
    parser.add_argument('-s', help='Overlapped MGEs list', dest='overlap_list', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    config['INPUT'] = dict()
    config['INPUT']['OVERLAP_LIST'] = args.overlap_list
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(args.reads_file)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(args.reads_file))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'],
                                                 config['INPUT']['INPUT_FILE_NAME_EXT'])
    config['OUTPUT'] = dict()
    config['OUTPUT']['OUT_DIR'] = os.path.abspath(args.output_dir_path)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

    # amr_set, mge_set, sample_colocalizations = get_colocalizations(config,
    #                                                                args.reads_file,
    #                                                                args.megares_sam,
    #                                                                args.mges_sam,
    #                                                                args.kegg_sam)
    amr_set, mge_set, sample_colocalizations = get_colocalizations(config,
                                                                   args.reads_file,
                                                                   args.megares_sam,
                                                                   args.mges_sam)

    logger = logging.getLogger()
    logger.info("Found {} reads with colocalizations".format(len(sample_colocalizations)))

    # Output colocalizations to stdout
    csvfile = sys.stdout
    writer = csv.writer(csvfile)
    # header = ['read', 'ARG', 'ARG positions', 'MGE(s)', 'MGE(s) positions', 'KEGG(s)', 'KEGG(s) positions']
    header = ['read', 'ARG', 'ARG positions', 'MGE(s)', 'MGE(s) positions']
    writer.writerow(header)
    for read, colocalization in sample_colocalizations.items():
        row = list()
        row.append(read)
        # ARG
        arg = ""
        arg_position = ""
        for gene in colocalization:
            if gene[2] == 'amr':
                if arg != "":
                    arg += ';'
                    arg_position += ';'
                arg += gene[0]
                arg_position += "{}:{}".format(gene[1][0], gene[1][1])

        # Mges
        mges = ""
        mges_positions = ""
        for gene in colocalization:
            if gene[2] == 'mge':
                if mges != "":
                    mges += ';'
                    mges_positions += ';'
                mges += gene[0]
                mges_positions += "{}:{}".format(gene[1][0], gene[1][1])

        # Kegg
        # keggs = ""
        # keggs_positions = ""
        # for gene in colocalization:
        #     if gene[2] == 'kegg':
        #         if keggs != "":
        #             keggs += ';'
        #             keggs_positions += ';'
        #         keggs += gene[0]
        #         keggs_positions += "{}:{}".format(gene[1][0], gene[1][1])

        # row.extend([arg, arg_position, mges, mges_positions, keggs, keggs_positions])
        row.extend([arg, arg_position, mges, mges_positions])
        writer.writerow(row)


if __name__ == "__main__":
    main()
