
# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

from Bio import SearchIO

from common import *


def main():
    parser = argparse.ArgumentParser(description='Cluster reads based on read length')
    parser.add_argument('-p', help='BLAT alignment file (.pls)', type=str, dest='pls_file', required=True)
    parser.add_argument('-s', help='Similarity threshold', type=float, dest='similarity_threshold', required=True)
    args = parser.parse_args()

    root_logger = init_logger()

    # find duplicates
    duplication_sets = []  # container of sets, where each set represents a clique of CCS that closely align
    # each alignment in psl: if qualifications met, those two belong to the same set
    cur_pls_file = args.pls_file
    root_logger.info('Reading in ' + cur_pls_file)
    for qresult in SearchIO.parse(cur_pls_file, 'blat-psl'):
        for hit in qresult:
            for hsp in hit.hsps:
                if sum(hsp.hit_span_all) >= args.similarity_threshold * hit.seq_len and sum(
                        hsp.query_span_all) >= args.similarity_threshold * qresult.seq_len:
                    if not duplication_sets:
                        duplication_sets.append(set())
                        duplication_sets[0].add(qresult.id)
                        duplication_sets[0].add(hit.id)
                    else:
                        alignment_is_in_set = False
                        for dup_set in duplication_sets:
                            if (qresult.id in dup_set):
                                dup_set.add(hit.id)
                                alignment_is_in_set = True
                                break
                            elif (hit.id in dup_set):
                                dup_set.add(qresult.id)
                                alignment_is_in_set = True
                                break
                        if not alignment_is_in_set:
                            duplication_sets.append(set())
                            duplication_sets[-1].add(qresult.id)
                            duplication_sets[-1].add(hit.id)

    non_singleton_dup_sets = [dup_set for dup_set in duplication_sets if len(dup_set) > 1]
    # append dup set to tsv
    tsv_handle = sys.stdout
    tsv_writer = csv.writer(tsv_handle)
    tsv_writer.writerows(sorted(non_singleton_dup_sets, key=lambda dup_set: len(dup_set), reverse=True))

if __name__ == "__main__":
    main()
