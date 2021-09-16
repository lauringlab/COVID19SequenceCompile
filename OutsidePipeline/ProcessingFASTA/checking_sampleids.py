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
meta_file_location = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20210908_Nanopore_Run_48/"
meta_file = "20210908_Nanopore_Run_48.meta.csv"
meta = pd.read_csv(meta_file_location + meta_file, index_col = None, header = 0, dtype = object)

plate_sampleids = meta['sample_id'].tolist()

## second, get a list of all the sample_ids in the fasta file we're comparing it to
fasta_file = "20210908_Nanopore_Run_48.all.consensus.renamed.full.fasta"
file_1 = meta_file_location + fasta_file

fasta_file_sampleids = list()

for record in SeqIO.parse(file_1, "fasta"):

    id = str(record.id).split()[0]
    id = id.split("/", 1)[0] ## removes excess info that comes through with barcode in some instances (NB01/ARTIC/nanopolish)

    fasta_file_sampleids.append(id)


## print out -- all matches OR some different (list them) OR all different
