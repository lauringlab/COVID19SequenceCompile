

### Author: Andrew Valesano
### Project: SARS-CoV-2 Sequencing
### Purpose: Change the names of the fasta entries for uploading to GISAID.

# Usage: Run in the batch folder in ProcessedGenomes.

#"""
# use in git bash
#c:/users/juliegil/appdata/local/programs/python/python38/python.exe C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/ProcessingFASTA/prep_fasta_for_gisaid.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20210517_Nanopore_Run_25/20210517_Nanopore_Run_25"


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

    if ("juliegil" in os.getcwd()):
        s_path = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/"
        s_path2 = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_B/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/"
        #full_loc = "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("leighbak" in os.getcwd()):
        s_path = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/"
        #full_loc = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    else:
        print("Current working directory username not recognized.")

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    #meta_file = s_path + args.prefix + ".forgisaid.meta.csv"

    # read in meta file, which is the compiled file (full_compiled_data.csv)
    #meta = pd.read_csv(meta_file, index_col = None, header = 0, dtype = object)
    #meta.columns = meta.columns.astype(str)

    file_1 = s_path + "gisaid_rsv.fasta"
    file_2 = s_path + "gisaid_back_info.csv"
    file_2b = s_path2 + "gisaid_back_info.csv"
    #file_3 = s_path + args.prefix + ".all.consensus.final.gisaid.fasta"

    all_fasta = list()
    for record in SeqIO.parse(file_1, "fasta"):

        id = str(record.id)
        all_fasta.append(id)

    fasta_headers = pd.DataFrame (all_fasta, columns = ['headers'])

    # split headers on document
    # words = string.split(',')
    # user_df['name'].str.split(pat = ' ', expand = True)
    # >hRSV/A/USA/MA-IVY-D13J37N6/2022|EPI_ISL_17367997|2022-11-04

    fasta_headers[['organism', 'subtype', 'country', 'id_string', 'year_tag_date']] = fasta_headers['headers'].str.split(pat = "/", expand = True)
    fasta_headers[['year', 'epi_isl', 'date']] = fasta_headers['year_tag_date'].str.split(pat = "|", expand = True)
    fasta_headers[['state', 'project', 'sample_id']] = fasta_headers['id_string'].str.split(pat = "-", expand = True)

    # Write
    #with open(file_2, 'w') as corrected:
    #    SeqIO.write(all_fasta, corrected, "fasta")
    fasta_headers.to_csv(file_2, encoding='utf-8', index=False)
    fasta_headers.to_csv(file_2b, encoding='utf-8', index=False)


if __name__ == "__main__":
    main()
