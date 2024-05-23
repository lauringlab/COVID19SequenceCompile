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

# ========================= Functions =============================

# ========================= Main =============================

def main():

    if ("juliegil" in os.getcwd()):
        s_path = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/"
        full_loc = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("leighbak" in os.getcwd()):
        s_path = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/"
        full_loc = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"
    elif ("chbl" in os.getcwd()):
        s_path = "/Users/chbl/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/"
        full_loc = "/Users/chbl/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"    
    else:
        print("Current working directory username not recognized.")

    parser = argparse.ArgumentParser()
    # takes the string after "--prefix" in the command line for use later
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    pn = args.prefix
    pn_date = pn[0:8]

    sequence_folder = s_path + pn + "/Segment_sequences/"
    onlyfiles = [f for f in os.listdir(sequence_folder) if os.path.isfile(os.path.join(sequence_folder, f))]

    #print(onlyfiles[1])

    files_avail = pd.DataFrame(onlyfiles, columns = ['file_name'])
    files_avail[['segment_info', 'fasta_ending']] = files_avail['file_name'].str.split('\\.', expand=True)
    files_avail[['sample_id', 'typeAorB', 'segment']] = files_avail['segment_info'].str.split('_', expand=True)

    files_avail['count'] = files_avail.groupby('sample_id')['sample_id'].transform('count')
    # remove anything that doesn't have 8 segments
    files_avail = files_avail[files_avail["count"] == 8]
    files_avail.fillna("",inplace=True)
    files_avail['true_segment'] = files_avail["segment"] #+ files_avail["segment2"]
    #print(files_avail.head())
    # now we have a dataframe with all the available segment file information, split so we can make new header names.
    # we need to cross reference that against our metadata file list, so we know which ones we want to submit.

    #

    ### This section converts barcode (NBXX) into sample_id.

    # create set of file names
    #file_1 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.final.fasta"
    #file_2 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.tmp.fasta"
    #file_3 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.final.genbank.fasta" # This is the file to use in genbank
    meta_file = s_path + args.prefix + "/" + args.prefix + ".forgenbank.meta.csv" # This is made in genbank file prep code - just uploaded names + sample ifs

    # read in .meta.csv file, which is the compiled file (full_compiled_data.csv)
    meta = pd.read_csv(meta_file, index_col = None, header = 0, dtype = object)
    meta.columns = meta.columns.astype(str)

    # merge files_avail and meta, only keeping entries that are in files_avail
    to_keep = files_avail.merge(meta, how='inner', on='sample_id')
    to_keep['final_VirusName'] = to_keep['VirusName'] + to_keep['true_segment']
    print(to_keep.head())


    # Change NBXX into sample_id
    all_fasta = list()
    # read in .all.consensus.fasta, which is the output from the lab (fastq -> fasta)
    # that contains all the sequences from a plate run
    # replace sampleid as ID part with the corresponding Virus name that is submitted to gisaid
    #print(range(0, to_keep.shape[0]))
    for i in range(0, to_keep.shape[0]):
        file_1 = s_path + args.prefix + "/Segment_sequences/" + to_keep['file_name'].iloc[i]
        #print(file_1)
        new_id_name = str(to_keep['final_VirusName'].iloc[i])
        print(new_id_name)
        for record in SeqIO.parse(file_1, "fasta"):

            record.id = ""
            record.description = str(new_id_name) + " [BioProject=PRJNA1023536]"
            #print(record.id)
            all_fasta.append(record)

    #print(len(all_fasta))
    # Write everything out as .all.consensus.tmp.fasta, with the replaced system
    file_2 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.final.genbanktmp.fasta"
    file_3 = s_path + args.prefix + "/" + args.prefix + ".all.consensus.final.genbank.fasta"
    #print(file_2)
    with open(file_3, 'w') as corrected:
        SeqIO.write(all_fasta, corrected, "fasta")

    #sed_cmd = """ sed '/^>/ s/ .*//' """ + '"' + file_2 + '"' + " > " + '"' + file_3 + '"'
    #os.system(sed_cmd)
    ## select all the IDs that we care about:
    #IDs = meta['sample_id'].tolist()
    #print(IDs)
    # Change NBXX into sample_id
    #all_fasta = list()
    # read in .all.consensus.fasta, which is the output from the lab (fastq -> fasta)
    # that contains all the sequences from a plate run
    # replace sampleid as ID part with the corresponding Virus name that is submitted to gisaid
    #for record in SeqIO.parse(file_1, "fasta"):

    #    id = str(record.id).split()[0]
    #    id = id.split("/", 1)[0] ## removes excess info that comes through with barcode in some instances (NB01/ARTIC/nanopolish)
        #id = str(id).upper()

    #    if(id in IDs):
    #        meta_sample = meta[meta.sample_id == id]
    #        new_ID = list(set(meta_sample.VirusName))[0]

    #        record.id = str(new_ID) # needed to change numeric sample_ids to be recognized as character strings
    #        all_fasta.append(record)

    # Write everything out as .all.consensus.tmp.fasta, with the replaced system
    #with open(file_2, 'w') as corrected:
    #    SeqIO.write(all_fasta, corrected, "fasta")
    # rename to .all.consensus.final.gisaid.fasta
    # and remove everything after a space character on fasta entry lines. The Biopython modules add extra characters after the sample ID
    #sed_cmd = """ sed '/^>/ s/ .*//' """ + '"' + file_2 + '"' + " > " + '"' + file_3 + '"' # need quotes around file names with spaces in them (from the dropbox folder)
    #print(sed_cmd)
    #os.system(sed_cmd)
    #os.system("rm " + '"' + file_2 + '"') # remove the .all.consensus.tmp.fasta

if __name__ == "__main__":
    main()
