

### Author: Andrew Valesano
### Project: SARS-CoV-2 Sequencing
### Purpose: Change the names of the fasta entries for uploading to GISAID.

# Usage: Run in the batch folder in ProcessedGenomes.

"""

 python prep_fasta_for_gisaid.py --prefix 20210519_Nanopore_Run_26

"""

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# ========================= Functions =============================

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    ### This section converts barcode (NBXX) into sample_id.

    file_1 = args.prefix + ".all.consensus.fasta"
    file_2 = args.prefix + ".all.consensus.tmp.fasta"
    file_3 = args.prefix + ".all.consensus.renamed.full.fasta" # This is the file to use in pangolin and nextclade.
    meta_file = args.prefix + ".meta.csv" # This is the full compiled list, or the subset that matches this run.

    meta = pd.read_csv(meta_file, index_col = None, header = 0, dtype = object)
    meta.columns = meta.columns.astype(str)

    # Change NBXX into sample_id
    all_fasta = list()
    for record in SeqIO.parse(file_1, "fasta"):

        id = str(record.id).split()[0]
        print(id)

        meta_sample = meta[meta.SampleBarcode == id]
        new_ID = list(set(meta_sample.sample_id))[0]

        record.id = new_ID
        all_fasta.append(record)

    # Write everything
    with open(file_2, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")
    sed_cmd = """ sed '/^>/ s/ .*//' """ + file_2 + " > " + file_3
    print(sed_cmd)
    os.system(sed_cmd)
    os.system("rm " + file_2)

    ### Make filtered fasta file ready for GISAID. Length filtering step occurred in Nanopore script in lab (90% completeness, or 27000 bases).

    file_1 = args.prefix + ".all.consensus.final.fasta"
    file_2 = args.prefix + ".all.consensus.final.tmp.fasta"
    file_3 = args.prefix + ".all.consensus.final.gisaid.fasta"

    all_fasta = list()
    for record in SeqIO.parse(file_1, "fasta"):

        id = str(record.id).split()[0]
        print(id)

        meta_sample = meta[meta.sample_id == id]
        new_ID = list(set(meta_sample.strain))[0]

        record.id = new_ID
        all_fasta.append(record)

    # Write
    with open(file_2, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")
    sed_cmd = """ sed '/^>/ s/ .*//' """ + file_2 + " > " + file_3
    print(sed_cmd)
    os.system(sed_cmd)
    os.system("rm " + file_2)

if __name__ == "__main__":
    main()
