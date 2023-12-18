# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

import json
from io import StringIO
import argparse
import configparser
import time
import ast
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.KEGG import REST
from tqdm import tqdm


# Functions
def to_df(result):
    return pd.read_table(StringIO(result), header=None)

def download_database(organism_list, out_file_path):
    with open(out_file_path, 'w') as output_handle:
        for organism in organism_list:
            print(organism)
            genes_list_res = REST.kegg_list(organism).read()
            genes_list = list(to_df(genes_list_res)[0])
            for gene in tqdm(genes_list):
                fasta_seq = ""
                for wait_time in [0.00001, 16, 32, 64, 128, 256, 512]:
                    try:
                        time.sleep(wait_time)
                        fasta_seq = REST.kegg_get(gene, "ntseq").read()
                    except:
                        print("Error in downloading gene:\t{}\twating:{}".format(gene, wait_time))
                        continue
                    break
                if len(fasta_seq) == 0:
                    print("Error persists, exiting")
                    exit()
                with StringIO(fasta_seq) as fasta_io:
                    records = list(SeqIO.parse(fasta_io, "fasta"))
                    record = SeqRecord(Seq(str(records[0].seq).upper()), id=records[0].id, description='')
                    SeqIO.write(record, output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description='Cluster reads based on read length')
    parser.add_argument('-o', help='Output file.', type=str, dest='out_file', required=True)
    parser.add_argument('-g', help='KEGG organisms list.', nargs='+', dest='organisms_list', required=True)
    args = parser.parse_args()

    # Organisms List
    print(args.organisms_list)

    # Download KEGG's prokaryotes default database
    download_database(args.organisms_list, args.out_file)


if __name__ == "__main__":
    main()