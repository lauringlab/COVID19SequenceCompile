# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def main():


    # get every folder name
    root='/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes'
    dirlist = [ item for item in os.listdir(root) if os.path.isdir(os.path.join(root, item)) ]
    #dirlist.sort()
    #dirlist = dirlist[len(dirlist)-184:len(dirlist)]

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
        writeout_path = "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate/complete_set.all.consensus.renamed.full.fasta"
        with open(writeout_path, 'w') as corrected:
            SeqIO.write(all_fasta, corrected, "fasta")




if __name__ == "__main__":
    main()
