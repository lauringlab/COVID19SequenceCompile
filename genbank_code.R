################################################################################
#          Creation of GenBank Dataset for COVID-19 Genetic Sampling           #
#                         Last Updated: 08/30/2023                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
#library(openxlsx)

################################################################################
#                   Genbank Files - Upload and Data Checks                     #
################################################################################

# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/genbank")

### output location of genbank files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.tsv", path = genbank_fp)

genbank_storage <- data.frame()

for (i in file_list){
    genbank_in <- read.delim(paste0(genbank_fp, "/", i), header = TRUE, row.names = NULL)
    genbank_in <- genbank_in[, c(1:3)]
    colnames(genbank_in) <- c("Accession", "SequenceID", "Release")
    genbank_storage <- rbind(genbank_storage, genbank_in)
}


genbank_storage <- genbank_storage %>% distinct()  

# create sample_id column
genbank_storage$sample_id <-  sapply(strsplit(as.character(genbank_storage$SequenceID),'/'), "[", 1)

genbank_storage <- separate(genbank_storage, col=sample_id, into=c('one', 'two'), sep='-')

genbank_storage$sample_id <- genbank_storage$two

genbank_storage <- genbank_storage %>% select(SequenceID, Accession, sample_id)

### rename columns 
rename_columns <- c("genbank_SequenceID", "genbank_Accession", "sample_id")
colnames(genbank_storage) <- rename_columns

# write out the compiled file
write.csv(genbank_storage, paste0(outputLOC, "/sample_full_genbank_list.csv"), row.names = FALSE, na = "")
