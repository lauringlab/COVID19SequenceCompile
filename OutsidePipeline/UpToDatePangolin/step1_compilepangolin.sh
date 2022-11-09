#!/bin/bash

UNIQNAME="$1"

if [ "$UNIQNAME" = "juliegil" ]; then
  python3 /Users/$UNIQNAME/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/OutsidePipeline/UpToDatePangolin/compile_all_fastas.py
else
  python3 /Users/$UNIQNAME/Documents/Lauring_Lab/COVID19SequenceCompile/OutsidePipeline/UpToDatePangolin/compile_all_fastas.p
fi

cd '/Users/$UNIQNAME/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate'
scp complete_set.all.consensus.renamed.full.fasta $UNIQNAME@greatlakes-xfer.arc-ts.umich.edu:/nfs/turbo/umms-alauring/pangolin_refreshes/full_fasta_file_folder
