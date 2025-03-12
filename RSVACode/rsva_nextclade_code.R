################################################################################
#        Creation of NextClade Dataset for RSVA Genetic Sampling           #
#                         Last Updated: 01/17/2025                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
#library(openxlsx)

################################################################################
#                 NextClade Files - Upload and Data Checks                     #
################################################################################

# nextclade file path
nc_fp <- paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/SequenceOutcomes/nextclade")

### output location of nextclade files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list_rsva <- list.files(pattern = "*.tsv", path = nc_fp)

nc_storage <- data.frame()
nc_old <- data.frame()
nc_new <- data.frame()

for (each_file in file_list_rsva){
  if(grepl("2023", each_file) | grepl("2024", each_file)){
    nc1 <- read.table(paste0(nc_fp, "/", each_file), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
    nc1 <- nc1 %>% select(seqName, clade, G_clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
    colnames(nc1) <- c("seqName", "clade", "G_clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs")
    nc_old <- rbind(nc1, nc_old)
  } else {
    nc2 <- read.table(paste0(nc_fp, "/", each_file), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
    nc2 <- nc2 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
    nc2$G_clade <- ""
    nc2 <- nc2 %>% select(seqName, clade, G_clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
    colnames(nc2) <- c("seqName", "clade", "G_clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs")
    nc_new <- rbind(nc2, nc_new)
  }
  
  nc_storage <- rbind(nc_old, nc_new)
  
  nc_storage$nextclade_runDate <- date_from_file_FIRST(each_file)
}
### remove commas from clade column
nc_storage$clade <- gsub(",", "", nc_storage$clade)

### rename columns 
rename_columns <- c("SampleID", "nextclade_clade", "nextclade_Gclade", "nextclade_totalMissing", "nextclade_qcOverallScore", "nextclade_qcOverallStatus", "nextclade_totalMutations", "nextclade_totalNonACGTNs", "nextclade_runDate")
colnames(nc_storage) <- rename_columns

################################################################################
## calculate genome completeness. genome size = 15,225 (gisaid reference sequence)
# changing to 15225 - 400 = 14,825
# Completeness = 100*(genome_size - totalMissing)/genome_size

nc_storage$nextclade_completeness <- 100*(14825 - as.numeric(nc_storage$nextclade_totalMissing)) / 14825

################################################################################

# write out the compiled file
write.csv(nc_storage, paste0(outputLOC, "/sample_full_nextclade_list.csv"), row.names = FALSE, na = "")
