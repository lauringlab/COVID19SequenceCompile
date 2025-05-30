################################################################################
#        Creation of NextClade Dataset for hMPV Genetic Sampling           #
#                         Last Updated: 01/17/2025                           #
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
nc_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/SequenceOutcomes/nextclade")

### output location of nextclade files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.tsv", path = nc_fp)

nc_storage <- data.frame()

for (each_page in file_list){
  nc1 <- read.table(paste0(nc_fp, "/", each_page), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
  
  #print(colnames(nc1))
  if("totalMutations" %in% colnames(nc1)){
    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalMutations, totalNonACGTNs, aaSubstitutions)
  } else {
    ### nextclade update changed column totalMutations to totalSubstitutions (near 6/18/2021)
    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs, aaSubstitutions)
    colnames(nc1) <- c("seqName", "clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs", "aaSubstitutions")
  }
  
  
  ### add date column from file name
  nc1$nextclade_runDate <- date_from_file_FIRST(each_page)
  
  nc_storage <- rbind(nc_storage, nc1)
}

### remove commas from clade column
nc_storage$clade <- gsub(",", "", nc_storage$clade)

### rename columns 
rename_columns <- c("SampleID", "nextclade_clade", "nextclade_totalMissing", "nextclade_qcOverallScore", "nextclade_qcOverallStatus", "nextclade_totalMutations", "nextclade_totalNonACGTNs", "SF456L_present", "nextclade_runDate")
colnames(nc_storage) <- rename_columns

################################################################################
## calculate genome completeness. genome size = 13350 (NCBI ref seq length)
# Completeness = 100*(genome_size - totalMissing)/genome_size

nc_storage$nextclade_completeness <- 100*(13350 - as.numeric(nc_storage$nextclade_totalMissing)) / 13350

################################################################################

# write out the compiled file
write.csv(nc_storage, paste0(outputLOC, "/sample_full_nextclade_list.csv"), row.names = FALSE, na = "")
