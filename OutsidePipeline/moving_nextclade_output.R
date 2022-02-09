library(tidyverse)
library(openxlsx)

options(scipen = 999)

# get nextclade.tsv from nextclade run
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2
# \3_ProcessedGenomes\{plate name}
# and drop it in
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\
# 4_SequenceSampleMetadata\SequenceOutcomes\nextclade
# renamed

if (file.exists(paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/nextclade.tsv"))){
      file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/nextclade.tsv"),
                to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/nextclade/", plate_name, "_nextclade.tsv"), 
                overwrite = TRUE)
} else {

      file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/", plate_name, "_nextclade.tsv"),
                to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/nextclade/", plate_name, "_nextclade.tsv"), 
                overwrite = TRUE)
}