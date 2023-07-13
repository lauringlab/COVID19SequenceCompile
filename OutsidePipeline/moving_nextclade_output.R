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

if (grepl("SC2", plate_name)){
    if (file.exists(paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/nextclade.tsv"))){
          file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/nextclade.tsv"),
                    to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/nextclade/", plate_name, "_nextclade.tsv"), 
                    overwrite = TRUE)
    } else {
    
          file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/", plate_name, "_nextclade.tsv"),
                    to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/nextclade/", plate_name, "_nextclade.tsv"), 
                    overwrite = TRUE)
    }
  
} else if (grepl("IAV", plate_name)){
  # check for H1 file
  if (file.exists(paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_name, "/", plate_name, "_HA_H1_nextclade.tsv"))){
    
    file.copy(from = paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_name, "/", plate_name, "_HA_H1_nextclade.tsv"),
              to = paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/nextclade/", plate_name, "_HA_H1_nextclade.tsv"), 
              overwrite = TRUE)
    
    
  } 
  
  # check for H3 file
  if (file.exists(paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_name, "/", plate_name, "_HA_H3_nextclade.tsv"))){
    
    file.copy(from = paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_name, "/", plate_name, "_HA_H3_nextclade.tsv"),
              to = paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/nextclade/", plate_name, "_HA_H3_nextclade.tsv"), 
              overwrite = TRUE)
  } 
  
}

