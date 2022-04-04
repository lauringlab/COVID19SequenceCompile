################################################################################
#          Creation of GISAID Dataset for COVID-19 Genetic Sampling            #
#                         Last Updated: 05/18/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
#library(openxlsx)

################################################################################
#                   GISAID Files - Upload and Data Checks                      #
################################################################################

# gisaid file path
gisaid_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/gisaid")

### output location of gisaid files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.tsv", path = gisaid_fp)

## check here - there should only be one .tsv file in this folder location, the most
# current gisaid file
if (length(file_list) != 1){
  stop(paste0("There is more than one .tsv file in ", gisaid_fp))
}

gisaid_storage <- read.delim(paste0(gisaid_fp, "/", file_list[1]))
  
# remove any empty rows/columns that may come in
gisaid_storage <- remove_empty(gisaid_storage)
  
# select columns we care about
gisaid_storage <- gisaid_storage %>% select(Virus.name, Accession.ID)

# create sample_id column
gisaid_storage$sample_id <-  sapply(strsplit(as.character(gisaid_storage$Virus.name),'/'), "[", 3)
gisaid_storage <- filter(gisaid_storage, grepl("MI-UM", Virus.name) | grepl("IVY", Virus.name) | grepl("RVTN", Virus.name))
gisaid_storage$sample_id <-  sapply(strsplit(as.character(gisaid_storage$sample_id),'-'), "[", 3)

### rename columns 
rename_columns <- c("gisaid_strain", "gisaid_epi_isl", "sample_id")
colnames(gisaid_storage) <- rename_columns

# write out the compiled file
write.csv(gisaid_storage, paste0(outputLOC, "/sample_full_gisaid_list.csv"), row.names = FALSE, na = "")
