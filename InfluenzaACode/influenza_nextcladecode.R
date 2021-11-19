################################################################################
#      Creation of NextClade Dataset for INFLUENZA A Genetic Sampling          #
#                           Created: 11/19/2021                                #
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
nc_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/nextclade")

### output location of nextclade files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

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
  
  nc1$nextclade_type <- trimws(as.character(strsplit(each_page, "_")[[1]][7]))
  
  nc_storage <- rbind(nc_storage, nc1)
}


### rename columns 
rename_columns <- c("SampleID", "nextclade_HA_clade", "nextclade_HA_totalMissing", "nextclade_HA_qcOverallScore", "nextclade_HA_qcOverallStatus", "nextclade_HA_totalMutations", "nextclade_HA_totalNonACGTNs", "nextclade_HA_runDate", "nextclade_HA_type")
colnames(nc_storage) <- rename_columns

################################################################################
## calculate genome completeness. 
# Completeness = 100*(genome_size - totalMissing)/genome_size

#summary(as.numeric(nc_storage$nextclade_HA_totalMissing))

# H3N2 = 1737
# H1N1 = 1752
# if only coding region, both = 1701

#nc_storage$nextclade_completeness <- 100*(29903 - as.numeric(nc_storage$nextclade_totalMissing)) / 29903

################################################################################

# write out the compiled file
write.csv(nc_storage, paste0(outputLOC, "/sample_full_nextclade_list.csv"), row.names = FALSE, na = "")
