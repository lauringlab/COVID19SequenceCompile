#### Get list of manifest files
#### purpose: determine which manifest files are new relative to last run
#### and decrease processing time

#### Julie Gilbert
#### 2022-10-07


library(tidyverse)
library(lubridate)

# set manifest file location
manifest_folder_location <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests")

# get a list of all folders that contain manifests
all_manifest_folders <- list.dirs(manifest_folder_location, recursive = FALSE)
# recursive must be false so that the folders within each folder are not listed
# we don't want to go into those in this process

# remove the ManifestsComplete entry
all_manifest_folders <- all_manifest_folders[all_manifest_folders != "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/ManifestsComplete"]

# iterate through each folder
# and generate a list of all file names inside

all_manifest_files_list <- c()

for (each_folder in all_manifest_folders){
  # get all .csv files
  files_inside <- list.files(each_folder, include.dirs = FALSE, recursive = FALSE, full.names = FALSE, pattern = "*.csv")
  all_manifest_files_list <- c(all_manifest_files_list, files_inside)
  # get all .xlsx files
  files_inside <- list.files(each_folder, include.dirs = FALSE, recursive = FALSE, full.names = FALSE, pattern = "*.xlsx")
  all_manifest_files_list <- c(all_manifest_files_list, files_inside)
}

#### write out the list of file names
saveRDS(all_manifest_files_list, paste0(manifest_folder_location, "/ManifestsComplete/current_manifest_list.RDS"))




