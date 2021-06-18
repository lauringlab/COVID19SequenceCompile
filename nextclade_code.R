################################################################################
#        Creation of NextClade Dataset for COVID-19 Genetic Sampling           #
#                         Last Updated: 05/18/2021                             #
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
nc_fp <- paste0(starting_path, "SequenceSampleMetadata/SequenceOutcomes/nextclade")

### output location of nextclade files, all together
outputLOC <- paste0(starting_path, "SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.tsv", path = nc_fp)

nc_storage <- data.frame()

for (each_page in file_list){
  nc1 <- read.table(paste0(nc_fp, "/", each_page), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
  
  #print(colnames(nc1))
  if("totalMutations" %in% colnames(nc1)){
    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalMutations, totalNonACGTNs)
  } else {
    ### nextclade update changed column totalMutations to totalSubstitutions (near 6/18/2021)
    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
    colnames(nc1) <- c("seqName", "clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs")
  }
  
  ### add date column from file name
  nc1$nextclade_runDate <- date_from_file_FIRST(each_page)
  
  nc_storage <- rbind(nc_storage, nc1)
}


### rename columns 
rename_columns <- c("SampleID", "nextclade_clade", "nextclade_totalMissing", "nextclade_qcOverallScore", "nextclade_qcOverallStatus", "nextclade_totalMutations", "nextclade_totalNonACGTNs", "nextclade_runDate")
colnames(nc_storage) <- rename_columns

################################################################################
## calculate genome completeness. genome size = 29903
# Completeness = 100*(genome_size - totalMissing)/genome_size

nc_storage$nextclade_completeness <- 100*(29903 - as.numeric(nc_storage$nextclade_totalMissing)) / 29903

################################################################################

# write out the compiled file
write.csv(nc_storage, paste0(outputLOC, "/sample_full_nextclade_list.csv"), row.names = FALSE, na = "")
