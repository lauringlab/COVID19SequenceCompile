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
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_name, "/")

# read in compiled dataset
seq_list <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################
# filter to plate run
seq_list2 <- filter(seq_list, PlateNumber == strsplit(plate_name, "_")[[1]][5] & PlateDate == paste0(substr(strsplit(plate_name, "_")[[1]][1], 1, 4), "-", substr(strsplit(plate_name, "_")[[1]][1], 5, 6), "-", substr(strsplit(plate_name, "_")[[1]][1], 7, 8)))
#puis <- filter(seq_list, grepl("pui", tolower(flag)) | grepl("pui", tolower(SampleSourceLocation)))

missing_subject_id <- nrow(filter(seq_list2, subject_id == "" | is.na(subject_id)))
missing_collection_date <- nrow(filter(seq_list2, coll_date == "" | is.na(coll_date)))
missing_received_source <- nrow(filter(seq_list2, received_source == "" | is.na(received_source)))

if (missing_subject_id > 3){
  print("Check manifests for:")
  print(plate_name)
  stop("Subject IDs Missing")
}

if (missing_collection_date > 3){
  print("Check manifests for:")
  print(plate_name)
  stop("Collection Dates Missing")
}

if (missing_received_source > 3){
  print("Check manifests for:")
  print(plate_name)
  stop("Received Sources Missing")
}

if (nrow(seq_list2) != 96){
  print("There are not 96 samples on this plate:")
  print(plate_name)
  print(paste0("Plate rows = ", nrow(seq_list2)))
}

# for use for identifying missing manifests
# out <- filter(seq_list2, subject_id == "")
# write.csv(out, "C:/Users/juliegil/Dropbox (University of Michigan)/Personal_DropBox/2021/MissingManifests/run40_20210823.csv", row.names = FALSE, na = "")

# write out that file as the .meta.csv file - change name as appropriate
write.csv(seq_list2, paste0(outputLOC, plate_name, ".meta.csv"), row.names = FALSE, na = "")
