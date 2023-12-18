# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

from common import *

from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Cluster reads based on read length')
    parser.add_argument('-r', help='Reads File (FASTA/FASTQ)', type=str, dest='reads_file', required=True)
    parser.add_argument('-d', help='Csv with duplicates sets', type=str, dest='duplicates_file', required=True)
    args = parser.parse_args()

    root_logger = init_logger()

    fastq_with_dups = args.reads_file
    sets_csv = args.duplicates_file

    sets = []
    with open(sets_csv, 'r') as sets_handle:
        sets_csv_reader = csv.reader(sets_handle)
        for row in sets_csv_reader:
            curr_set_elements = []
            for read_id in row:
                curr_set_elements.append(read_id)
            sets.append(set(curr_set_elements))

    # Sets are not disjoint :(
    are_sets_disjoint = False
    itr = 0
    while (not are_sets_disjoint):
        # print("iteration: " + str(itr))
        # print(len(sets))
        itr += 1
        are_sets_disjoint = True
        combined_sets = []
        used_set_indices = {}
        for i in range(0, len(sets)):
            used_set_indices[i] = False

        for i in range(0, len(sets)):
            if not used_set_indices[i]:
                combined_set = sets[i]
                used_set_indices[i] = True
            else:
                # print("continue to next from " + str(i))
                continue

            for j in range(i + 1, len(sets)):
                if not used_set_indices[j]:
                    if not len(combined_set & sets[j]) == 0:
                        # print("merging " + str(i) + ", " + str(j))
                        combined_set = combined_set | sets[j]
                        are_sets_disjoint = False
                        used_set_indices[j] = True
                else:
                    continue

            # print("appending: " + str(i))
            combined_sets.append(combined_set)

        # print(len(combined_sets))
        sets = combined_sets

    set_deduplicated = {}
    frozen_sets = []
    # print(len(sets))
    for curr_set in sets:
        frozen_set = frozenset(curr_set)
        frozen_sets.append(frozen_set)
        set_deduplicated[frozen_set] = False

    if (is_gz_file(fastq_with_dups)):
        root_logger.info('Opening gzipped file')
        file_handler = gzip.open(fastq_with_dups, 'rt')
    else:
        root_logger.info('Opening uncompressed file')
        file_handler = open(fastq_with_dups, 'rt')

    dedup_records = []
    set_sizes = []
    num_singletons = 0
    for record in SeqIO.parse(file_handler, "fastq"):
        record_in_set = False
        for dup_set in frozen_sets:
            if record.id in dup_set:
                record_in_set = True
                curr_record_set = dup_set
                break

        if record_in_set:
            if not set_deduplicated[curr_record_set]:
                dedup_records.append(record)
                set_sizes.append(len(curr_record_set))
                set_deduplicated[curr_record_set] = True
        else:  # Singletons are not in the set csv
            num_singletons += 1
            dedup_records.append(record)

    file_handler.close()
    out_handle = sys.stdout
    SeqIO.write(dedup_records, out_handle, "fastq")

if __name__ == "__main__":
    main()