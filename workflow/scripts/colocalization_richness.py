#!/usr/bin/env python

# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

from common import *

#Essentially just a way to define equality/uniqueness
class Colocalization:
    distance_cutoff = 250
    separator = ":::"

    def __init__(self, amr_group ="", mge_header ="", distance = -1):
        self.amr_group = amr_group
        self.mge_header = mge_header
        self.distance = distance

    def __str__(self):
        return self.amr_group + Colocalization.separator + self.mge_header + Colocalization.separator + str(self.distance)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.amr_group == other.amr_group and \
               self.mge_header == other.mge_header and \
               abs(self.distance - other.distance) <= Colocalization.distance_cutoff

    def __lt__(self, other):
        if self.amr_group < other.amr_group:
            return True
        elif self.mge_header < other.mge_header:
            return True
        else:
            return self.distance < other.distance

    def __le__(self, other):
        if self.amr_group <= other.amr_group:
            return True
        elif self.mge_header <= other.mge_header:
            return True
        else:
            return self.distance <= other.distance


    def __gt__(self, other):
        if self.amr_group > other.amr_group:
            return True
        elif self.mge_header > other.mge_header:
            return True
        else:
            return self.distance > other.distance


    def __ge__(self, other):
        if self.amr_group >= other.amr_group:
            return True
        elif self.mge_header >= other.mge_header:
            return True
        else:
            return self.distance >= other.distance



def main():
    parser = argparse.ArgumentParser(description='Colocalizations Finder.')
    parser.add_argument('-i', help='Colocalizations csv file', dest='coloc_csv_path', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    root_logger = init_logger()

    #Create ontology dictionary from MEGARes ontology file
    megares_ontology, _ = read_megares_v2_ontology(config)

    colocs = []
    coloc_counts = {}
    with open(args.coloc_csv_path, 'r') as colocalizations_csv:
        colocalizations_reader = csv.reader(colocalizations_csv)
        header = next(colocalizations_reader)
        for row in colocalizations_reader:
            read_name = row[0]
            ARGs_list = row[1].split(";")
            ARGs_positions = row[2].split(";")
            MGEs_list = row[3].split(";")
            MGEs_positions = row[4].split(";")

            # Only rows with coloc
            if (len(ARGs_list) > 0 and len(MGEs_list) > 0) and (ARGs_list[0] != '' and MGEs_list[0] != ''):
                AMR_gene = ARGs_list[0]
                AMR_gene_start = int(ARGs_positions[0].split(":")[0])
                AMR_gene_end = int(ARGs_positions[0].split(":")[1])
                AMR_group = megares_ontology[AMR_gene]["group"]

                for idx, MGE_gene in enumerate(MGEs_list):
                    MGE_gene_start = int(MGEs_positions[idx].split(":")[0])
                    MGE_gene_end = int(MGEs_positions[idx].split(":")[1])

                    if AMR_gene_end < MGE_gene_start:
                        distance = int(MGE_gene_start) - int(AMR_gene_end)
                    else :
                        distance = int(AMR_gene_start) - int(MGE_gene_end)

                    coloc = Colocalization(AMR_group, MGE_gene, distance)

                    ref_coloc = Colocalization()
                    is_new_coloc = True
                    for unique_coloc in colocs:
                        if unique_coloc == coloc:
                            ref_coloc = unique_coloc
                            is_new_coloc = False
                            break

                    if is_new_coloc:
                        colocs.append(coloc)
                        coloc_counts[coloc] = 1
                    else:
                        coloc_counts[ref_coloc] += 1

    sorted_unique_colocs = sorted(colocs, key=lambda coloc: coloc_counts[coloc], reverse=True)

    root_logger.info("Found {} unique colocalizations".format(len(sorted_unique_colocs)))

    richness_csv = sys.stdout
    coloc_writer = csv.writer(richness_csv)
    coloc_writer.writerow([])
    coloc_writer.writerow(["Number of unique colocalizations:", len(sorted_unique_colocs)])
    coloc_writer.writerow([])
    coloc_writer.writerow(["MEGARes group", "MGE gene", "Distance (within " + str(Colocalization.distance_cutoff) + " nts)", "Occurences"])
    for coloc in sorted_unique_colocs:
        coloc_writer.writerow([coloc.amr_group, coloc.mge_header, coloc.distance, coloc_counts[coloc]])


if __name__ == "__main__":
    main()

