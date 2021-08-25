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
outputLOC <- paste0(starting_path, "ProcessedGenomes/20210816_Nanopore_Run_40/")

# read in compiled dataset
seq_list <- read.csv(paste0(starting_path, "SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################
# filter to plate run
seq_list2 <- filter(seq_list, PlateNumber == "40" & PlateDate == "2021-08-16")

# for use for identifying missing manifests
# out <- filter(seq_list2, subject_id == "")
# write.csv(out, "C:/Users/juliegil/Dropbox (University of Michigan)/Personal_DropBox/2021/MissingManifests/run40_20210823.csv", row.names = FALSE, na = "")

# write out that file as the .meta.csv file - change name as appropriate
write.csv(seq_list2, paste0(outputLOC, "20210818_Nanopore_Run_41.meta.csv"), row.names = FALSE, na = "")
