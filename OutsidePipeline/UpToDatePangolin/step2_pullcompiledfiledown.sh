#!/bin/bash

UNIQNAME="$1"
TODAYDATE=$(date "+%Y%m%d")

#echo $TODAYDATE
# move into dropbox folder where complete fasta information is
cd '/Users/'$UNIQNAME'/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate'

# move the old full pangolin file into archive
mv lineage_report*.csv '/Users/'$UNIQNAME'/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate/ARCHIVE/'

# copy a new full pangolin file back down
scp $UNIQNAME@greatlakes-xfer.arc-ts.umich.edu:/nfs/turbo/umms-alauring/pangolin_refreshes/pan_out/lineage_report.csv lineage_report_$TODAYDATE.csv
