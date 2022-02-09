### Author: Andrew Valesano/Julie (Jules) Gilbert
### Project: SARS-CoV-2 Sequencing
### Purpose: Change the names of the fasta entries for uploading to GISAID

# Usage: Run in the batch folder in ProcessedGenomes.

#"""
# use in git bash
#c:/users/juliegil/appdata/local/programs/python/python38/python.exe C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/ProcessingFASTA/prep_fasta_NumberTwo.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20220131_SC2_Nanopore_Run_110/20220131_SC2_Nanopore_Run_110"
# python prep_fasta_NumberTwo.py --prefix 20210519_Nanopore_Run_26
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
    file_1 = args.prefix + ".all.consensus.renamed.full.fasta"
    file_2 = args.prefix + ".all.consensus.tmp.fasta"
    file_3 = args.prefix + ".all.consensus.final.gisaid.fasta" # This is the file to use in gisaid
    meta_file = args.prefix + ".forgisaid.meta.csv" # This is made in gisaid file prep code - just uploaded names + sample ifs

    # read in .meta.csv file, which is the compiled file (full_compiled_data.csv)
    meta = pd.read_csv(meta_file, index_col = None, header = 0, dtype = object)
    meta.columns = meta.columns.astype(str)
    ## select all the IDs that we care about:
    IDs = meta['sample_id'].tolist()

    # Change NBXX into sample_id
    all_fasta = list()
    # read in .all.consensus.fasta, which is the output from the lab (fastq -> fasta)
    # that contains all the sequences from a plate run
    # replace sampleid as ID part with the corresponding Virus name that is submitted to gisaid
    for record in SeqIO.parse(file_1, "fasta"):

        id = str(record.id).split()[0]
        id = id.split("/", 1)[0] ## removes excess info that comes through with barcode in some instances (NB01/ARTIC/nanopolish)
        #id = str(id).upper()
        print(id)

        if(id in IDs):
            meta_sample = meta[meta.sample_id == id]
            new_ID = list(set(meta_sample.VirusName))[0]

            record.id = str(new_ID) # needed to change numeric sample_ids to be recognized as character strings
            all_fasta.append(record)

    # Write everything out as .all.consensus.tmp.fasta, with the replaced system
    with open(file_2, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")
    # rename to .all.consensus.final.gisaid.fasta
    # and remove everything after a space character on fasta entry lines. The Biopython modules add extra characters after the sample ID
    sed_cmd = """ sed '/^>/ s/ .*//' """ + '"' + file_2 + '"' + " > " + '"' + file_3 + '"' # need quotes around file names with spaces in them (from the dropbox folder)
    print(sed_cmd)
    os.system(sed_cmd)
    os.system("rm " + '"' + file_2 + '"') # remove the .all.consensus.tmp.fasta

if __name__ == "__main__":
    main()
