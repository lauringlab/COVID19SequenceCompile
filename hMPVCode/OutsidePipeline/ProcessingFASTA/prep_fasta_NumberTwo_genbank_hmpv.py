### Author: Andrew Valesano/Julie (Jules) Gilbert/Leigh Papalambros
### Project: hMPV Sequencing
### Purpose: Change the names of the fasta entries for uploading to Genbank after being processed through VADR

# Usage: Run in the batch folder in ProcessedGenomes.

#"""
# use in git bash
#c:/users/juliegil/appdata/local/programs/python/python38/python.exe C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/ProcessingFASTA/prep_fasta_NumberTwo.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20220208_SC2_Nanopore_Run_117/20220208_SC2_Nanopore_Run_117"
# python prep_fasta_NumberTwo.py --prefix "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20221222_SC2_Illumina_Run_75/20221222_SC2_Illumina_Run_75"
# python prep_fasta_NumberTwo.py --prefix 20210519_Nanopore_Run_26 --gb 0210519_nanopore_run_26
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

    if ("leighbak" in os.getcwd()):
        s_path = "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/3_ProcessedGenomes/"
        g_path = "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/6_GenBank_Uploads/"
        full_loc = "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("chbla" in os.getcwd()):
        s_path = "C:/Users/chbla/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/hMPV/3_ProcessedGenomes/"
        full_loc = "C:/Users/chbla/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/hPMV/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("chbl" in os.getcwd()):
        s_path = "/Users/chbl/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/3_ProcessedGenomes/"
        full_loc = "/Users/chbl/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"    
    else:
        print("Current working directory username not recognized.")

    parser = argparse.ArgumentParser()
    # takes the string after "--prefix" in the command line for use later
    parser.add_argument('--prefix', action="store", dest="prefix")
    #parser.add_argument('--gb', action="store", dest="gb")
    args = parser.parse_args()

    ### This section converts barcode (NBXX) into sample_id.

    # create set of file names
    file_1 = s_path + args.prefix + "/" + args.prefix + ".consensus.all.fasta"
    #file_2 = s_path + args.prefix + ".all.consensus.final.genbanktmp.fasta"
    file_3 = s_path + args.prefix + ".all.consensus.vadr.fasta" # This is the file to use in genbank
    file_4 = g_path + "upload_genbank_" + args.prefix + "/" + args.prefix + ".all.consensus.vadr.fasta" # This is the file to use in genbank
    meta_file = s_path + args.prefix + "/" + args.prefix + ".forvadr.meta.csv" # This is made in gisaid file prep code - just uploaded names + sample ifs

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
        #id = id.split("/", 1)[0] ## removes excess info that comes through with barcode in some instances (NB01/ARTIC/nanopolish)
        #id = str(id).upper()
        print(id)

        if(id in IDs):
            meta_sample = meta[meta.sample_id == id]
            new_ID = list(set(meta_sample.VirusName))[0]
            record.id = ""
            record.description = str(new_ID) + " [BioProject=PRJNA1181811]" # needed to change numeric sample_ids to be recognized as character strings
            all_fasta.append(record)

    # Write everything out as .all.consensus.tmp.fasta, with the replaced system
    file_2 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.vadr.temp.fasta"
    file_3 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.vadr.fasta"
    file_4 = g_path + "upload_genbank_"  + args.prefix + "/" + args.prefix + ".all.consensus.vadr.fasta"
    with open(file_3, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")
    with open(file_4, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")    
     #rename to .all.consensus.final.gisaid.fasta
    # and remove everything after a space character on fasta entry lines. The Biopython modules add extra characters after the sample ID
    #sed_cmd = """ sed '/^>/ s/ .*//' """ + '"' + file_2 + '"' + " > " + '"' + file_3 + '"' # need quotes around file names with spaces in them (from the dropbox folder)
    #print(sed_cmd)
    #os.system(sed_cmd)
    #os.system("rm " + '"' + file_2 + '"') # remove the .all.consensus.tmp.fasta

if __name__ == "__main__":
    main()
