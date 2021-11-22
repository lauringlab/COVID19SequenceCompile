## Checking fasta file matches to plates
## expectation = the same
### Jules

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Process = samples have barcode map that assigns sample_ids to each sample
## sequence. If the barcode map is incorrect (contains the wrong samples) then
## the sequences will be matched with the wrong samples

## first, get a list of all the sample_ids in the plate map we're checking
meta_file_location = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20211115_SC2_Nanopore_Run_74/"
meta_file = "20211115_SC2_Nanopore_Run_74.meta.csv"

# read in file as pandas table
meta = pd.read_csv(meta_file_location + meta_file, index_col = None, header = 0, dtype = object)

# pull out only the sample id column, and convert it to a list
plate_sampleids = meta['sample_id'].tolist()

## second, get a list of all the sample_ids in the fasta file we're comparing it to
fasta_file = "20211115_SC2_Nanopore_Run_74.all.consensus.renamed.full.fasta"
file_1 = meta_file_location + fasta_file

fasta_file_sampleids = list()

# parse through each record of the fasta file
for record in SeqIO.parse(file_1, "fasta"):

    id = str(record.id).split()[0]
    id = id.split("/", 1)[0] ## removes excess info that comes through with barcode in some instances (NB01/ARTIC/nanopolish)

    fasta_file_sampleids.append(id)

## get any differences
## print out -- all matches OR some different (list them) OR all different
print("Sample IDs on the Plate Map that are NOT in the FASTA file: ")
print(list(set(plate_sampleids) - set(fasta_file_sampleids)))

print("\n")

print("Sample IDs in the FASTA file that are NOT on the Plate Map: ")
print(list(set(fasta_file_sampleids) - set(plate_sampleids)))
