# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def main():

    current_user = os.getcwd()

    if "juliegil" in current_user:
        root='/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes'
        writeout_path = "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate/complete_set.all.consensus.renamed.full.fasta"
    elif "leighbaker" in current_user:
        root='/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes'
        writeout_path = "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate/complete_set.all.consensus.renamed.full.fasta"
    elif "leighbak" in current_user:
        root='/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes'
        writeout_path = "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate/complete_set.all.consensus.renamed.full.fasta"
    else:
        print("User not recognized.")

    # get every folder name

    dirlist = [ item for item in os.listdir(root) if os.path.isdir(os.path.join(root, item)) ]
    #dirlist.sort()
    #dirlist = dirlist[len(dirlist)-184:len(dirlist)]
    dirlist.remove('20220720_SC2_Illumina_Run_58')
    dirlist.remove('20220722_SC2_Illumina_Run_60')
    dirlist.remove('20220817_SC2_Illumina_Run_63')
    dirlist.remove('20220823_SC2_Illumina_Run_64')
    dirlist.remove('20220926_SC2_Illumina_Run_67')
    dirlist.remove('20221004_SC2_Illumina_Run_68')
    dirlist.remove('20230124_SC2_Illumina_Run_78')
    dirlist.remove('20230330_SC2_Illumina_Run_95')
    dirlist.remove('20230523_SC2_Illumina_Run_106')
    #print(len(dirlist))
    all_fasta = list()
    for each_folder in dirlist:
        x = 0
        complete_path = root + "/" + each_folder + "/" + each_folder + ".all.consensus.renamed.full.fasta"

        try:
            for record in SeqIO.parse(complete_path, "fasta"):
                all_fasta.append(record)
            x = 1
        except:
            print("No File for " + each_folder)

    if x == 1:

        with open(writeout_path, 'w') as corrected:
            SeqIO.write(all_fasta, corrected, "fasta")




if __name__ == "__main__":
    main()
