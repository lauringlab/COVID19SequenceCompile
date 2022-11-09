#!/bin/bash

UNIQNAME="$1"

# run code to create full compiled fasta file
if [ "$UNIQNAME" = "juliegil" ]; then
  python3 /Users/$UNIQNAME/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/OutsidePipeline/UpToDatePangolin/compile_all_fastas.py
else
  python3 /Users/$UNIQNAME/Documents/Lauring_Lab/COVID19SequenceCompile/OutsidePipeline/UpToDatePangolin/compile_all_fastas.py
fi

# move into dropbox location where full fasta file is written to
cd '/Users/$UNIQNAME/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate'
# copy that full fasta file up into greatlakes
scp complete_set.all.consensus.renamed.full.fasta $UNIQNAME@greatlakes-xfer.arc-ts.umich.edu:/nfs/turbo/umms-alauring/pangolin_refreshes/full_fasta_file_folder
