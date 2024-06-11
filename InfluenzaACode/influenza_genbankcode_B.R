################################################################################
#          Creation of GenBank Dataset for Influenza Genetic Sampling          #
#                         Last Updated: 02/07/2024                            #
#                 Code Edited By: Leigh Papalambros                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
#library(openxlsx)

################################################################################
#                   Genbank Files - Upload and Data Checks                     #
################################################################################

# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/SequenceOutcomes/genbank")

### output location of genbank files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.tsv", path = genbank_fp)

genbank_storage <- data.frame()

for (i in file_list){
  genbank_in <- read.delim(paste0(genbank_fp, "/", i), header = TRUE, row.names = NULL)
  
  # get submission id
  a <- strsplit(i, "\\_")[[1]][2]
  b <- strsplit(a, "\\.")[[1]][1]
  
  genbank_in <- genbank_in[, c(1:3)]
  colnames(genbank_in) <- c("Accession", "SequenceID", "Release")
  
  # need to edit sequenceID to separate sample id out from the segment type
  genbank_in <- genbank_in %>% separate(SequenceID, into = c("id"), sep = "-", remove = FALSE)
  
  genbank_in$interest <- substring(genbank_in$SequenceID, regexpr("-", genbank_in$SequenceID) + 1, nchar(genbank_in$SequenceID))
  
  genbank_in <- genbank_in %>% mutate(segment_type = case_when(grepl("MP", interest) ~ "genbank_MP",
                                                               grepl("NP", interest) ~ "genbank_NP", 
                                                               grepl("NS", interest) ~ "genbank_NS", 
                                                               grepl("PA", interest) ~ "genbank_PA",
                                                               grepl("PB1", interest) ~ "genbank_PB1",
                                                               grepl("PB2", interest) ~ "genbank_PB2",
                                                               grepl("HA", interest) ~ "genbank_HA", 
                                                               grepl("NA", interest) ~ "genbank_NA",
                                                               T ~ "unknown"), 
                                      sample_id = case_when(grepl("MP", interest) ~ gsub("MP", "", interest),
                                                            grepl("NP", interest) ~ gsub("NP", "", interest), 
                                                            grepl("NS", interest) ~ gsub("NS", "", interest), 
                                                            grepl("PA", interest) ~ gsub("PA", "", interest),
                                                            grepl("PB1", interest) ~ gsub("PB1", "", interest),
                                                            grepl("PB2", interest) ~ gsub("PB2", "", interest),
                                                            grepl("HA", interest) ~ gsub("HA", "", interest),
                                                            grepl("NA", interest) ~ gsub("NA", "", interest),
                                                            T ~ "unknown"))
  

  genbank_in$SubmissionID <- b
  
  genbank_storage <- rbind(genbank_storage, genbank_in)
}


genbank_storage <- genbank_storage %>% distinct()  

####################

genbank_storage <- genbank_storage %>% select(id, Accession, sample_id, segment_type, SubmissionID)

genbank_storage <- reshape2::dcast(genbank_storage, SubmissionID + sample_id + id ~ segment_type, value.var = c("Accession"))

# # create sample_id column
# genbank_storage$sample_id <-  sapply(strsplit(as.character(genbank_storage$SequenceID),'/'), "[", 1)
# 
# genbank_storage <- separate(genbank_storage, col=sample_id, into=c('one', 'two'), sep='-')
# 
# genbank_storage <- genbank_storage %>% mutate(loc_code2 = one)
# 
# genbank_storage$sample_id <- genbank_storage$two

#genbank_storage <- genbank_storage %>% select(SequenceID, SubmissionID, Accession, sample_id, loc_code2)

### rename columns 
rename_columns <- c("genbank_SubmissionID", "sample_id", "loc_code2", "genbank_HA", "genbank_MP",
                    "genbank_NA", "genbank_NP", "genbank_NS", "genbank_PA", "genbank_PB1",
                    "genbank_PB2")
colnames(genbank_storage) <- rename_columns



# write out the compiled file
write.csv(genbank_storage, paste0(outputLOC, "/sample_full_genbank_list.csv"), row.names = FALSE, na = "")

