# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

import json
import os.path

from Bio import SeqIO
from sklearn.cluster import KMeans

from common import *


def main():
    parser = argparse.ArgumentParser(description='Cluster reads based on read length')
    parser.add_argument('-r', help='Reads File (FASTA/FASTQ)', type=str, dest='reads_file', required=True)
    parser.add_argument('-o', help='Output directory', type=str, dest='out_dir', required=True)
    parser.add_argument('-n', help='Number of clusters to generate', type=int, dest='num_clusters', default=200)
    parser.add_argument('-l', help='Reads Length File', type=str, dest='reads_lengths', required=True)
    args = parser.parse_args()

    root_logger = init_logger()

    reads_file = args.reads_file
    output_dir = args.out_dir
    num_clusters = args.num_clusters

    with open(args.reads_lengths) as json_file:
        reads_lengths = json.load(json_file)

    # Cluster read lengths
    reads_lengths_list = ([], [])
    for read_name, read_length in reads_lengths.items():
        reads_lengths_list[0].append(read_length)
        reads_lengths_list[1].append(read_name)

    X = np.array(reads_lengths_list[0])

    reads_clusters = dict()
    if len(reads_lengths.values()) < num_clusters:
        root_logger.error("Number of clusters selected ({}) is more than the reads in the input ({}).".format(num_clusters, len(reads_lengths.values())))
        root_logger.error("All reads in first cluster.")
        for read_name, read_length in reads_lengths.items():
            reads_clusters[read_name] = 0
    else:
        kmeans = KMeans(n_clusters=num_clusters).fit(X.reshape((X.shape[0], 1)))
        # Associate each read name to its cluster
        reads_clusters = dict()
        for j in range(0, kmeans.labels_.shape[0]):
            reads_clusters[reads_lengths_list[1][j]] = kmeans.labels_[j]

    # Open clusters file
    clusters_file_handles = list()
    for i in range(num_clusters):
        output_filename = output_dir + '/' + os.path.splitext(os.path.basename(reads_file))[0] + '_' + str(
            i) + '.fasta.gz'
        try:
            clusters_file_handles.append(gzip.open(output_filename, "wt"))
        except Exception as e:
            root_logger.error("Error while opening {}. Exception:".format(output_filename, e))
            exit(1)

    # Iterate through every input read, write to appropriate cluster
    if is_gz_file(reads_file):
        root_logger.info('Opening gzipped file')
        file_handler = gzip.open(reads_file, 'rt')
    else:
        root_logger.info('Opening uncompressed file')
        file_handler = open(reads_file, 'rt')

    for record in SeqIO.parse(file_handler, "fastq"):
        if (record.name in reads_clusters):
            cluster_label = reads_clusters[record.name]
        else:
            root_logger.error("Read {} not found in read lenghts file".format(record.names))
            exit(1)

        SeqIO.write(record, clusters_file_handles[cluster_label], 'fasta')

    for i in range(num_clusters):
        clusters_file_handles[i].close()


if __name__ == "__main__":
    main()
