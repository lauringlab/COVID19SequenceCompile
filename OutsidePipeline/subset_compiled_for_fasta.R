################################################################################
#       Creation of Subsetting Compiled File for FASTA Name Replacement        #
#                         Last Updated: 06/02/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

# load libraries
library(tidyverse)
library(lubridate)

################################################################################
# set paths 
starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
## set run folder accordingly
outputLOC <- paste0(starting_path, "ProcessedGenomes/20210803_Nanopore_Run_37/")

# read in compiled dataset
seq_list <- read.csv(paste0(starting_path, "SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################
# filter to plate run
seq_list2 <- filter(seq_list, PlateNumber == "37")
# write out that file as the .meta.csv file - change name as appropriate
write.csv(seq_list2, paste0(outputLOC, "20210803_Nanopore_Run_37.meta.csv"), row.names = FALSE, na = "")
