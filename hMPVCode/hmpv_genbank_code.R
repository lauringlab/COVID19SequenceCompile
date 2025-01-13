################################################################################
#          Creation of GenBank Dataset for RSVA Genetic Sampling               #
#                         Last Updated: 12/08/2023                             #
#                 Code Edited By: Leigh Papalambros                            #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
library(readr)
#library(openxlsx)

################################################################################
#                   Genbank Files - Upload and Data Checks                     #
################################################################################

# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/SequenceOutcomes/genbank")

### output location of genbank files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.txt", path = genbank_fp)

#table_test <- read.table(paste0(genbank_fp, "/seqids_20231130.txt"), header = TRUE)
#test_table_2 <- read.table(paste0(genbank_fp, "/seqids_20231208.txt"))

genbank_storage <- data.frame()

for (i in file_list){
  genbank_in <- read.table(paste0(genbank_fp, "/", i), header = TRUE)
  genbank_storage <- rbind(genbank_storage, genbank_in) 
}

#sapply(genbank_storage, class)



#for (i in file_list){
#  genbank_in <- read.table(paste0(genbank_fp, "/", i), header = TRUE, row.names = NULL)
  
  # get submission id
 # a <- strsplit(i, "\\_")[[1]][2]
#  b <- strsplit(a, "\\.")[[1]][1]
  
 # genbank_in <- genbank_in[, c(1:3)]
#  colnames(genbank_in) <- c("SubmissionID", "SequenceID", "Accession")
  
  #genbank_in$SubmissionID <- b
  
 # genbank_storage <- rbind(genbank_storage, genbank_in)
#}


genbank_storage <- genbank_storage %>% distinct()  

#sapply(genbank_storage, class)
# create sample_id column
genbank_storage$sample_id <- sapply(strsplit(as.character(genbank_storage$SampleID),'-'),"[", 3)


############
#genbank_storage <- separate(genbank_storage, col=SequenceID, into=c('one', 'two', "three"), sep='-')

#genbank_storage <- genbank_storage %>% mutate(loc_code2 = one)

#genbank_storage$sample_id <- genbank_storage$three

#genbank_storage <- genbank_storage %>% select(SequenceID, SubmissionID, Accession, sample_id, loc_code2)

### rename columns 
rename_columns <- c("genbank_SubmissionID", "genbank_SequenceID",  "genbank_Accession", "sample_id") #, "loc_code2")
colnames(genbank_storage) <- rename_columns



# write out the compiled file
write.csv(genbank_storage, paste0(outputLOC, "/sample_full_genbank_list.csv"), row.names = FALSE, na = "")




##########################
#testing section for parsing the genbank .txt files
