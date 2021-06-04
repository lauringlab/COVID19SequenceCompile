################################################################################
#        Creation of Pangolin Dataset for COVID-19 Genetic Sampling            #
#                         Last Updated: 05/19/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
#library(openxlsx)

################################################################################
#                 Pangolin Files - Upload and Data Checks                      #
################################################################################

# pangolin file path
pangolin_fp <- paste0(starting_path, "SequenceSampleMetadata/SequenceOutcomes/pangolin")

### output location of pangolin files, all together
outputLOC <- paste0(starting_path, "SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.csv", path = pangolin_fp)

pang_storage <- data.frame()

counter <- 1

for (each_page in file_list){
  pang1 <- read.csv(paste0(pangolin_fp, "/", each_page), colClasses = "character")
  #print(colnames(pang1))
  # remove any empty rows/columns that may come in
  pang1 <- remove_empty(pang1)

  pang1$pang1_runDate <- date_from_file_FIRST(each_page)
  
  
  if (counter == 1){
    pang_storage <- pang1
  } else {
    pang_storage <- merge(pang_storage, pang1, all = TRUE)
  }
  
  counter <- counter + 1
}

#colnames(pang_storage)

pang_storage <- pang_storage %>% select(taxon, lineage, probability, pangoLEARN_version, status, note, conflict, pango_version, pangolin_version, pang1_runDate) %>% distinct()
### rename columns 
rename_columns <- c("SampleID", "pangolin_lineage", "pangolin_probability", "pangoLEARN_version", "pangolin_status", "pangolin_note", "pangolin_conflict", "pango_version", "pangolin_version", "pangolin_runDate")
colnames(pang_storage) <- rename_columns

# write out the compiled file
write.csv(pang_storage, paste0(outputLOC, "/sample_full_pangolin_list.csv"), row.names = FALSE, na = "")
