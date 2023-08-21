

### Author: Andrew Valesano
### Project: SARS-CoV-2 Sequencing
### Purpose: Change the names of the fasta entries for uploading to GISAID.

# Usage: Run in the batch folder in ProcessedGenomes.

#"""
# use in git bash
#c:/users/juliegil/appdata/local/programs/python/python38/python.exe C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/ProcessingFASTA/prep_fasta_for_gisaid.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20210517_Nanopore_Run_25/20210517_Nanopore_Run_25"


# python3 prep_fasta_for_gisaid_rsva.py --prefix "20230323_RSVA_Illumina_Run_3/20230323_RSVA_Illumina_Run_3"

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

    if ("juliegil" in os.getcwd()):
        #s_path = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/3_ProcessedGenomes/"
        s_path = "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/3_ProcessedGenomes/"
        #full_loc = "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("leighbak" in os.getcwd()):
        s_path = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/3_ProcessedGenomes/"
        #full_loc = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    else:
        print("Current working directory username not recognized.")

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    meta_file = s_path + args.prefix + "/" + args.prefix + ".forgenbank.meta.csv"
    #meta_file = "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/3_ProcessedGenomes/20230307_RSVA_Illumina_Run_1/20230307_RSVA_Illumina_Run_1.forgenbank.meta.csv"

    print(meta_file)
    # read in meta file, which is the compiled file (full_compiled_data.csv)
    meta = pd.read_csv(meta_file, index_col = None, header = 0, dtype = object)
    meta.columns = meta.columns.astype(str)

    #file_1 = s_path + args.prefix + ".90.consensus.fasta"
    file_1 = s_path + args.prefix + "/" + args.prefix + ".full.consensus.fasta"
    file_2 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.final.tmp.fasta"
    file_3 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.final.genbank.fasta"

    all_fasta = list()
    for record in SeqIO.parse(file_1, "fasta"):

        id = str(record.id).split()[0]
        #print(id)

        meta_sample = meta[meta.sample_id == id]
        if not meta_sample.empty:

            #meta_sample = meta[meta.sample_id == id]
            new_ID = list(set(meta_sample.VirusName))[0]

            record.id = str(new_ID) # needed to change numeric sample_ids to be recognized as character strings
            print(record.id)
            all_fasta.append(record)

    # Write
    with open(file_2, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")

    sed_cmd = """ sed '/^>/ s/ .*//' """ + '"' + file_2 + '"' + " > " + '"' + file_3 + '"' # need quotes around file names with spaces in them (from the dropbox folder)
    print(sed_cmd)
    os.system(sed_cmd)
    os.system("rm " + '"' + file_2 + '"') # remove the .all.consensus.tmp.fasta

if __name__ == "__main__":
    main()
