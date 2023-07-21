################################################################################
#       Creation of Subsetting Compiled File for Checking Completion           #
#                            Updated: 2023-07-21                               #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

# load libraries
library(tidyverse)
library(lubridate)
library(gt)

if (grepl("IAV", plate_name)){
################################################################################
# set paths 
#starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
## set run folder accordingly
outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_name, "/")

# read in compiled dataset
seq_list <- read.csv(paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")
} else {
  ## set run folder accordingly
  outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/", plate_name, "/")
  
  # read in compiled dataset
  seq_list <- read.csv(paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")
}
################################################################################
# filter to plate run
seq_list2 <- filter(seq_list, PlateNumber == strsplit(plate_name, "_")[[1]][5] & PlateDate == paste0(substr(strsplit(plate_name, "_")[[1]][1], 1, 4), "-", substr(strsplit(plate_name, "_")[[1]][1], 5, 6), "-", substr(strsplit(plate_name, "_")[[1]][1], 7, 8)))
#puis <- filter(seq_list, grepl("pui", tolower(flag)) | grepl("pui", tolower(SampleSourceLocation)))

######################################
# check how current the samples are

sample_years <- unique(year(as_date(seq_list2$coll_date))) 
sample_years <- sample_years[!is.na(sample_years)] # have to remove NA from control well rows date listing
current_year <- year(Sys.Date())

if (any(!sample_years %in% current_year)){
  message("Samples with collection dates not in current year")
  show_samples <- filter(seq_list2, year(as_date(coll_date)) != current_year)
  show_samples <- show_samples %>% select(sample_id, subject_id, coll_date, received_source, received_date)
  print(show_samples %>% gt())
  letter_value <- readline(prompt="If you'd like to proceed, press y, if you'd like to stop, press n: ")
  if (letter_value[1] == "n"){
    stop("Stopped code to correct collection dates.")
  } 
}

######################################

missing_subject_id <- nrow(filter(seq_list2, subject_id == "" | is.na(subject_id)))
missing_collection_date <- nrow(filter(seq_list2, coll_date == "" | is.na(as_date(coll_date)) | is.na(coll_date)))
missing_received_source <- nrow(filter(seq_list2, received_source == "" | is.na(received_source)))


if (missing_subject_id > (96 - num_iavs_on_plate)){
  print("Check manifests for:")
  print(plate_name)
  stop("Subject IDs Missing")
}

if (missing_collection_date > (96 - num_iavs_on_plate)){
  print("Check manifests for:")
  print(plate_name)
  stop("Collection Dates Missing")
}

if (missing_received_source > (96 - num_iavs_on_plate)){
  print("Check manifests for:")
  print(plate_name)
  stop("Received Sources Missing")
}

if (nrow(seq_list2) != num_iavs_on_plate + 3){
  print(paste0("There are not ", num_iavs_on_plate, " samples on this plate:"))
  print(plate_name)
  print(paste0("Plate rows = ", nrow(seq_list2)))
}


# for use for identifying missing manifests
# out <- filter(seq_list2, subject_id == "")
# write.csv(out, "C:/Users/juliegil/Dropbox (University of Michigan)/Personal_DropBox/2021/MissingManifests/run40_20210823.csv", row.names = FALSE, na = "")

# write out that file as the .meta.csv file - change name as appropriate
write.csv(seq_list2, paste0(outputLOC, plate_name, ".meta.csv"), row.names = FALSE, na = "")
