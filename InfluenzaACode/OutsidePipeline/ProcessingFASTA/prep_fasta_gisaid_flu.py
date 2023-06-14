### Author: Andrew Valesano/Julie (Jules) Gilbert
### Project: SARS-CoV-2 Sequencing
### Purpose: Change the names of the fasta entries for uploading to Pangolin and Nextclade

# Usage: Run in the batch folder in ProcessedGenomes.

#"""
# use in git bash
#c:/users/juliegil/appdata/local/programs/python/python38/python.exe C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/ProcessingFASTA/prep_fasta_NumberOne.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/ProcessedGenomes/20210803_Nanopore_Run_37/20210803_Nanopore_Run_37"
# python prep_fasta_NumberOne.py --prefix 20210519_Nanopore_Run_26
# python prep_fasta_gisaid_flu.py --prefix 20230315_IAV_Illumina_Run_49
#"""

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import numpy as np
#from os import listdir
#from os.path import isfile, join

# ========================= Functions =============================

# ========================= Main =============================

def main():

    if ("juliegil" in os.getcwd()):
        s_path = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/"
        full_loc = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("leighbak" in os.getcwd()):
        s_path = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/"
        full_loc = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    else:
        print("Current working directory username not recognized.")

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    #pn_date = "20220926"
    #pn_run = "28"
    #pn = pn_date + "_IAV_Nanopore_Run_" + pn_run
    pn = args.prefix
    pn_date = pn[0:8]

    sequence_folder = s_path + pn + "/Segment_sequences/"
    onlyfiles = [f for f in os.listdir(sequence_folder) if os.path.isfile(os.path.join(sequence_folder, f))]

    loc90 = s_path + pn + "/"
    keep90s = list()
    #library_checker = {}
    file90 = loc90 + pn + ".90.consensus.fasta"
    for record in SeqIO.parse(file90, "fasta"):
        keep90s.append(record.id)
        #if record.seq not in library_checker:
            #library_checker[record.seq] = 1
        #else:
            #library_checker[record.seq] = library_checker[record.seq] + 1

    #for each_i in library_checker:
        #print(library_checker[each_i])

    df_out = pd.DataFrame(keep90s, columns=[""])
    df_out.to_csv(('{}gisaid_90keeps.csv').format(sequence_folder), index=False)

    seq_hide_rvtn = pd.read_csv(full_loc)
    seq_hide_rvtn = seq_hide_rvtn[['sample_id', 'sample_id_lauring']]
    # merge keeps90 version in df_out with seq_hide_rvtn
    sq_hide2 = pd.merge(df_out, seq_hide_rvtn, left_on = "", right_on = "sample_id")
    # if sample_id_lauring is not null, then keep that, otherwise keep the keep90s id
    sq_hide2['new_keeps'] = np.where(sq_hide2['sample_id_lauring'].notnull(), sq_hide2['sample_id_lauring'], sq_hide2[''])
    #print(sq_hide2.head())
    #keep90s = list(sq_hide2['new_keeps'])



    not_used = 0
    all_fasta = list()
    ids_only = list()
    for file_in in onlyfiles:
        file_next = sequence_folder + file_in
        #print(file_in)

        for record in SeqIO.parse(file_next, "fasta"):
            if str(record.id).split("_")[0] in keep90s:

                original = sq_hide2[sq_hide2['sample_id']==str(record.id).split("_")[0]]['new_keeps'].values[0] + "_" + record.id.split("_")[1]
                print(original)
                record.id = str(original) + "_" + pn_date
                all_fasta.append(">" + record.id)
                ids_only.append(record.id)
                all_fasta.append(record.seq)

            else:
                not_used = not_used + 1


    df_out = pd.DataFrame(all_fasta, columns=[""])

    df_out.to_csv(('{}gisaid_sequencelist.csv').format(sequence_folder), index=False)

    df_out2 = pd.DataFrame(ids_only, columns=["IDS"])
    df_out2.to_csv(('{}gisaid_IDlist.csv').format(sequence_folder), index=False)

    print("{} samples did not meet 90% rule.".format(not_used))


if __name__ == "__main__":
    main()
