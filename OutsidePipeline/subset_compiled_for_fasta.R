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
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/20211214_SC2_Nanopore_Run_88/")

# read in compiled dataset
seq_list <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################
# filter to plate run
seq_list2 <- filter(seq_list, PlateNumber == "88" & PlateDate == "2021-12-14")
#puis <- filter(seq_list, grepl("pui", tolower(flag)) | grepl("pui", tolower(SampleSourceLocation)))

# for use for identifying missing manifests
# out <- filter(seq_list2, subject_id == "")
# write.csv(out, "C:/Users/juliegil/Dropbox (University of Michigan)/Personal_DropBox/2021/MissingManifests/run40_20210823.csv", row.names = FALSE, na = "")

# write out that file as the .meta.csv file - change name as appropriate
write.csv(seq_list2, paste0(outputLOC, "20211214_SC2_Nanopore_Run_88.meta.csv"), row.names = FALSE, na = "")
