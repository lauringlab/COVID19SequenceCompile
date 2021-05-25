

### Author: Andrew Valesano
### Project: SARS-CoV-2 Sequencing
### Purpose: Select specific sequences by ID.

### Usage: python SelectSequences.py --infile in.fasta --outfile out.fasta --select IDs.csv
### IDs.csv must have a single column with selected fasta entries, with "ID" as the header.



# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', action="store", dest="infile")
    parser.add_argument('--select', action="store", dest="select")
    parser.add_argument('--outfile', action="store", dest="outfile")
    args = parser.parse_args()

    selected = pd.read_csv(args.select, index_col = None, header = 0, dtype = str)
    IDs = selected['ID'].tolist()

    all_fasta = list()
    for record in SeqIO.parse(args.infile, "fasta"):

        id = str(record.id).split()[0]
        #print(id)

        if(id in IDs):
            print(id)
            all_fasta.append(record)

    with open(args.outfile, 'w') as final:
        SeqIO.write(all_fasta, final, "fasta")

if __name__ == "__main__":
    main()
