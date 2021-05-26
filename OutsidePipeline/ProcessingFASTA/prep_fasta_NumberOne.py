### Author: Andrew Valesano/Julie (Jules) Gilbert
### Project: SARS-CoV-2 Sequencing
### Purpose: Change the names of the fasta entries for uploading to Pangolin and Nextclade

# Usage: Run in the batch folder in ProcessedGenomes.

#"""
# use in git bash
#c:/users/juliegil/appdata/local/programs/python/python38/python.exe C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/ProcessingFASTA/prep_fasta_for_gisaid.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/ProcessedGenomes/20210517_Nanopore_Run_25/20210517_Nanopore_Run_25"
# python prep_fasta_for_gisaid.py --prefix 20210519_Nanopore_Run_26
#"""

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
    # takes the string after "--prefix" in the command line for use later
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    ### This section converts barcode (NBXX) into sample_id.

    # create set of file names
    file_1 = args.prefix + ".all.consensus.fasta"
    file_2 = args.prefix + ".all.consensus.tmp.fasta"
    file_3 = args.prefix + ".all.consensus.renamed.full.fasta" # This is the file to use in pangolin and nextclade.
    meta_file = args.prefix + ".meta.csv" # This is the full compiled list, or the subset that matches this run.

    # read in .meta.csv file, which is the compiled file (full_compiled_data.csv)
    meta = pd.read_csv(meta_file, index_col = None, header = 0, dtype = object)
    meta.columns = meta.columns.astype(str)

    # Change NBXX into sample_id
    all_fasta = list()
    # read in .all.consensus.fasta, which is the output from the lab (fastq -> fasta)
    # that contains all the sequences from a plate run
    # replace barcode as ID part with the corresponding sample_id
    for record in SeqIO.parse(file_1, "fasta"):

        id = str(record.id).split()[0]
        print(id)

        meta_sample = meta[meta.SampleBarcode == id]
        new_ID = list(set(meta_sample.sample_id))[0]

        record.id = new_ID
        all_fasta.append(record)

    # Write everything out as .all.consensus.tmp.fasta, with the replaced system
    with open(file_2, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")
    # rename to .all.consensus.renamed.full.fasta
    sed_cmd = """ sed '/^>/ s/ .*//' """ + file_2 + " > " + file_3
    print(sed_cmd)
    os.system(sed_cmd)
    os.system("rm " + file_2) # remove the .all.consensus.tmp.fasta

if __name__ == "__main__":
    main()
