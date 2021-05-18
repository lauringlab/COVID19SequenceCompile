################################################################################
#        Creation of Pangolin Dataset for COVID-19 Genetic Sampling            #
#                         Last Updated: 05/18/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
#library(openxlsx)

################################################################################
#                 Pangolin Files - Upload and Data Checks                      #
################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

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
  
  # remove any empty rows/columns that may come in
  pang1 <- remove_empty(pang1)
  
  pang1_runDate <- trimws(as.character(strsplit(each_page, "_")[[1]][1]))
  pang1$pang1_runDate <- paste0(substr(pang1_runDate, 1, 4), "-", substr(pang1_runDate, 5, 6), "-", substr(pang1_runDate, 7, 8))
  
  
  if (counter == 1){
    pang_storage <- pang1
  } else {
    pang_storage <- merge(pang_storage, pang1, all = TRUE)
  }
  
  # check for column names: position, sample_id, subject_id, coll_date, flag
  #column_name_check <- colnames(pang1)
  #true_columns <- c("taxon", "lineage", "probability", "pangoLEARN_version", "status", "note")
  
  
  # if (identical(true_columns, column_name_check)){
  #   ## then do nothing
  #   x <- 0
  # } else {
  #   if (ncol(pang1) == 6){
  #     # if the number of columns is the same, just rename them
  #     colnames(pang1) <- true_columns
  #   } else {
  #     # find out which column is missing
  #     whatsdifferent <- setdiff(true_columns, column_name_check)
  #     print(whatsdifferent)
  #     print("There is a column difference in")
  #     print(each_page)
  #     stop()
  #   }
  # }
  
  #pang_storage <- rbind(pang_storage, pang1)
  
  counter <- counter + 1
}

#colnames(pang_storage)

pang_storage <- pang_storage %>% select(taxon, lineage, probability, pangoLEARN_version, status, note, conflict, pango_version, pangolin_version, pang1_runDate) %>% distinct()
### rename columns 
rename_columns <- c("SampleID", "pangolin_lineage", "pangolin_probability", "pangoLEARN_version", "pangolin_status", "pangolin_note", "pangolin_conflict", "pango_version", "pangolin_version", "pangolin_runDate")
colnames(pang_storage) <- rename_columns

# write out the compiled file
write.csv(pang_storage, paste0(outputLOC, "/sample_full_pangolin_list.csv"), row.names = FALSE, na = "")
