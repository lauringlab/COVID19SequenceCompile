library(tidyverse)
library(openxlsx)

options(scipen = 999)

# get lineage_report.csv from pangolin run
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2
# \3_ProcessedGenomes\{plate name}
# and drop it in
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\
# 4_SequenceSampleMetadata\SequenceOutcomes\pangolin
# renamed

file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/lineage_report.csv"),
          to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/", plate_name, "_pangolin.csv"), 
          overwrite = TRUE)
