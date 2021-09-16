## Checking fasta file matches to plates
## expectation = the same
import pandas as pd
## Process = samples have barcode map that assigns sample_ids to each sample
## sequence. If the barcode map is incorrect (contains the wrong samples) then
## the sequences will be matched with the wrong samples

## first, get a list of all the sample_ids in the plate map we're checking
meta_file_location = "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20210908_Nanopore_Run_48/"
meta_file = "20210908_Nanopore_Run_48.meta.csv"
meta = pd.read_csv(meta_file_location + meta_file, index_col = None, header = 0, dtype = object)

print(meta['sample_id'])
## second, get a list of all the sample_ids in the fasta file we're comparing it to

## print out -- all matches OR some different (list them) OR all different
