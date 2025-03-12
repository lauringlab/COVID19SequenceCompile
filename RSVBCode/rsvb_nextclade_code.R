################################################################################
#        Creation of NextClade Dataset for RSVB Genetic Sampling           #
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
nc_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/SequenceOutcomes/nextclade")

### output location of nextclade files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list_rsvb <- list.files(pattern = "*.tsv", path = nc_fp)

nc_storage <- data.frame()
nc_old <- data.frame()
nc_new <- data.frame()

for (each_file in file_list_rsvb){
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
#for (each_page in file_list){
#  nc1 <- read.table(paste0(nc_fp, "/", each_page), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
 
  #print(colnames(nc1))
#  if("G_clade" %in% colnames(nc1)){
#    nc1 <- nc1 %>% select(seqName, clade, G_clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
#    colnames(nc1) <- c("seqName", "clade", "G_clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs")
#  } else {
    ### nextclade update changed column totalMutations to totalSubstitutions (near 6/18/2021)
#    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
    #nc1_2$G_clade <- ""
    #nc1_2 <- nc1_2 %>% select(seqName, clade, G_clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs)
#    colnames(nc1) <- c("seqName", "clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs")
#  }
  #nc2 <- rbind(nc1_1, nc1_2)
  
  ### add date column from file name
#  nc1$nextclade_runDate <- date_from_file_FIRST(each_page)
  
 # nc_storage <- rbind(nc_storage, nc1)
#}




### remove commas from clade column
nc_storage$clade <- gsub(",", "", nc_storage$clade)

### rename columns 
rename_columns <- c("SampleID", "nextclade_clade", "nextclade_Gclade", "nextclade_totalMissing", "nextclade_qcOverallScore", "nextclade_qcOverallStatus", "nextclade_totalMutations", "nextclade_totalNonACGTNs", "nextclade_runDate")
colnames(nc_storage) <- rename_columns

################################################################################
## calculate genome completeness. genome size = 15,222 (gisaid reference sequence)
# Completeness = 100*(genome_size - totalMissing)/genome_size

nc_storage$nextclade_completeness <- 100*(15222 - as.numeric(nc_storage$nextclade_totalMissing)) / 15222

################################################################################

# write out the compiled file
write.csv(nc_storage, paste0(outputLOC, "/sample_full_nextclade_list.csv"), row.names = FALSE, na = "")

