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
gisaid_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/new_gisaid_try")

### output location of gisaid files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

################################################################################

file_list <- list.files(pattern = "*.tsv", path = gisaid_fp)

## check here - there should only be one .tsv file in this folder location, the most
# current gisaid file
#if (length(file_list) != 1){
#  stop(paste0("There is more than one .tsv file in ", gisaid_fp))
#}

gisaid_storage <- data.frame()

for (i in file_list){
    gisaid_in <- read.delim(paste0(gisaid_fp, "/", i))
    gisaid_storage <- rbind(gisaid_storage, gisaid_in)
}

# remove any empty rows/columns that may come in
#gisaid_storage <- remove_empty(gisaid_storage)
#only keeping distinct rows 10,000 row limit for Gisaid
gisaid_storage <- gisaid_storage %>% distinct()  
# select columns we care about
gisaid_storage <- gisaid_storage %>% select(Virus.name, Accession.ID, Clade, Lineage)

# create sample_id column
gisaid_storage$sample_id <-  sapply(strsplit(as.character(gisaid_storage$Virus.name),'/'), "[", 3)
#gisaid_storage <- filter(gisaid_storage, grepl("MI-UM", Virus.name) | grepl("IVY", Virus.name) | grepl("RVTN", Virus.name))
#gisaid_storage$sample_id2 <- sapply(strsplit(as.character(gisaid_storage$sample_id),'-'), "[", 3)
gisaid_storage <- separate(gisaid_storage, col=sample_id, into=c('one', 'two', 'three', 'four'), sep='-')
gisaid_storage$four[is.na(gisaid_storage$four)] <- ""
gisaid_storage <- gisaid_storage %>% mutate(sample_id = case_when(four == "" ~ paste0(gisaid_storage$three, gisaid_storage$four), 
                                                                  T ~ paste0(gisaid_storage$three, "-", gisaid_storage$four)))

gisaid_storage <- gisaid_storage %>% select(Virus.name, Accession.ID, Clade, Lineage, sample_id)

### rename columns 
rename_columns <- c("gisaid_strain", "gisaid_epi_isl", "gisaid_clade", "gisaid_pango_lineage", "sample_id")
colnames(gisaid_storage) <- rename_columns

# write out the compiled file
write.csv(gisaid_storage, paste0(outputLOC, "/sample_full_gisaid_list.csv"), row.names = FALSE, na = "")
