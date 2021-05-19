################################################################################
#         Creation of Manifest Dataset for COVID-19 Genetic Sampling           #
#                         Last Updated: 05/19/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
library(openxlsx)
library(withr)
source("pipeline_functions.R")

################################################################################
#                Manifest Files - Upload and Data Checks                       #
################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

# Manifest file paths (there should be a path per source)
cbr_manifest_fp <- paste0(starting_path, "SequenceSampleMetadata/Manifests/CBR")
martin_manifest_fp <- paste0(starting_path, "SequenceSampleMetadata/Manifests/Martin")
cstp_manifest_fp <- paste0(starting_path, "SequenceSampleMetadata/Manifests/CSTP")
edidnow_manifest_fp <- paste0(starting_path, "SequenceSampleMetadata/Manifests/EDIDNOW")
cdcivy_manifest_fp <- paste0(starting_path, "SequenceSampleMetadata/Manifests/CDCIVY")

manifest_folder_list <- c(cbr_manifest_fp, martin_manifest_fp, cstp_manifest_fp, 
                          edidnow_manifest_fp)

### output location of manifest files, all together
outputLOC <- paste0(starting_path, "SequenceSampleMetadata/Manifests/ManifestsComplete")

# each manifest provider will have their own folder, with all files inside
# don't want to re-build the entire thing every time, so will eventually need a 
# system that only pulls in the newest information, but the first pass of this 
# code will create the full set
# when "newest" added, create output text file report with # new items added
# any warnings (duplicates) should also print there

################################################################################

manifest_storage <- data.frame()
#duplicate_ssc <- data.frame()

# will iterate through folders
for (each_folder in manifest_folder_list){
    
    ### get names of all .csv files in folder
    file_list <- list.files(pattern = "*.csv", path = each_folder)
    
    # then iterate through files within each folder
    for (each_file in file_list){
      
        # read in the file
        file_in <- read.csv(paste0(each_folder, "/", each_file), colClasses = "character")
        
        # turn any "" or " " into NA
        file_in[file_in == ""] <- NA
        file_in[file_in == " "] <- NA
        
        # remove any empty rows that may come in
        file_in <- remove_empty(file_in, which = "rows")
        
        # check for column names: position, sample_id, subject_id, coll_date, flag
        column_name_check <- colnames(file_in)
        true_columns <- c("position", "sample_id", "subject_id", "coll_date", "flag")
        
        if (identical(true_columns, column_name_check)){
          ## then do nothing
          x <- 0
        } else {
          if (ncol(file_in) == 5){
            colnames(file_in) <- true_columns
          } else {
            # find out which column is missing
            whatsdifferent <- setdiff(true_columns, column_name_check)
            print(whatsdifferent)
            print("There is a column difference in")
            print(each_file)
            stop()
          }
        }
        
        # check that sample_id, subject_id, and coll_date are all filled in
        
        if(any(is.na(file_in$sample_id))){
          print(each_file)
          stop("There are missing sample ids.")
        }
        
        if(any(is.na(file_in$subject_id))){
          print(each_file)
          stop("There are missing subject ids.")
        }
        
        if(any(is.na(file_in$coll_date))){
          print(each_file)
          stop("There are missing collection dates.")
        }
        
        ## reformat coll_date to YYYY-MM-DD format if necessary
        test_date_format <- substr(as.character(file_in[1, 4]), 1, 4)
        
        if (is.na(as.numeric(test_date_format))){
             file_in$coll_date <- as.POSIXct(file_in$coll_date, format = "%m/%d/%y")
        }
        
        # add in 2 new columns: received_date and received_source (from file name)
        file_in$received_date <- date_from_file(each_file)
        
        rec_source <- trimws(as.character(strsplit(each_file, "_")[[1]][1]))
        file_in$received_source <- rec_source
        
        # bind all rows together
        manifest_storage <- rbind(manifest_storage, file_in)
        
    }
    
    ### select only distinct rows
    manifest_storage <- manifest_storage %>% distinct()
    
}


################################################################################
## handle cdc ivy manifests

# read in cdc ivy site code list 
cdc_sites <- read.csv(paste0(cdcivy_manifest_fp, "/Keys/CDC_SiteCodebook.csv"), colClasses = "character")

cdc_file_list <- list.files(pattern = "*.xlsx", path = cdcivy_manifest_fp)

cdc_ivy_storage <- data.frame()
full_ivy <- data.frame()

for (each_file in cdc_file_list){
  fileone <- read.xlsx(paste0(cdcivy_manifest_fp, "/", each_file), sheet = 1, detectDates = TRUE)
  fileone <- filter(fileone, !is.na(as.numeric(`Position.#`)))
  
  ## change all column names to lowercase and remove leading/lagging white space
  ## to make it easier to process
  cdc_names <- colnames(fileone)
  new_names <- c()
  for (i in cdc_names){
    new_names <- c(new_names, tolower(trimws(i)))
  }
  colnames(fileone) <- new_names
  
  ## keep full set of cdc ivy rows separate, for back checks on full data
  full_ivy <- rbind(full_ivy, fileone)
  
  fileone <- fileone %>% select(`position.#`, site.name, study.id, collection.date, aliquot.id)
  
  colnames(fileone) <- c("position", "SiteName", "subject_id", "coll_date", "sample_id")
  
  ### site name checks
  fileone$SiteName_check <- ifelse(fileone$SiteName %in% cdc_sites$SiteCode, 0, 1)
  
  if (sum(fileone$SiteName_check, na.rm = TRUE) != 0){
    print(each_file)
    stop("There are incorrect site names in the manifest.")
  } else {
    fileone <- fileone %>% select(position, SiteName, subject_id, coll_date, sample_id)
  }
  
  # add in 2 new columns: received_date and received_source (from file name)
  #rec_date <- trimws(as.character(strsplit(each_file, "_")[[1]][2]))
  #rec_date <- paste0(substr(rec_date, 1, 4), "-", substr(rec_date, 5, 6), "-", substr(rec_date, 7, 8))
  fileone$received_date <- date_from_file(each_file)
  
  rec_source <- trimws(as.character(strsplit(each_file, "_")[[1]][1]))
  fileone$received_source <- rec_source
  
  ### add in "regular" manifest columns
  fileone$flag <- NA
  
  ### re-arrange variables
  fileone <- fileone %>% select(position, sample_id, subject_id, coll_date, flag, received_date, received_source, SiteName)
  
  cdc_ivy_storage <- rbind(cdc_ivy_storage, fileone)
}

### write out full ivy set
write.csv(full_ivy, paste0(cdcivy_manifest_fp, "/Full_IVY_Set/IVY_sample_full_manifest_list.csv"), row.names = FALSE, na = "")

### add onto main manifest file HERE
manifest_storage$SiteName <- NA
manifest_storage <- rbind(manifest_storage, cdc_ivy_storage)


################################################################################
# check for sample_id/subject_id/coll_date duplicates
# count of unique sample_id, subject_id, coll_date combinations
unique_ids <- nrow(manifest_storage %>% select(sample_id, subject_id, coll_date) %>% distinct())

if (unique_ids != nrow(manifest_storage)){
  # identify duplicates
  dupes <- manifest_storage %>% group_by(sample_id, subject_id, coll_date) %>% summarize(count_unique = length(sample_id))
  # merge with the original file
  mfs <- merge(manifest_storage, dupes, by = c("sample_id", "subject_id", "coll_date"), all.x = TRUE)
  # filter that down to a set of duplicates, to use in the output report
  duplicate_ssc <- filter(mfs, count_unique != 1)
  ### alter flag column in original file to note the duplication
  manifest_storage$flag <- ifelse(is.na(manifest_storage$flag), "", manifest_storage$flag)
  manifest_storage$flag <- ifelse(manifest_storage$sample_id %in% duplicate_ssc$sample_id & manifest_storage$subject_id %in% duplicate_ssc$subject_id 
                                  & manifest_storage$coll_date %in% duplicate_ssc$coll_date, 
                                  paste0(manifest_storage$flag, " ", "Duplicate Sample - Subject - Collection"), manifest_storage$flag)
  manifest_storage$flag <- trimws(manifest_storage$flag)
  manifest_storage$flag <- ifelse(manifest_storage$flag == "", NA, manifest_storage$flag)
  
  #duplicate_ssc <- rbind(duplicate_ssc, duplicate_ssc2)
}


################################################################################
# check for instances where leading zeros may have been dropped from subject id
# CBR should == 9 digits (MRN)
# CSTP should == 8 digits (UMID)

# add a column to check length of subject_id
manifest_storage$subject_id_length <- nchar(manifest_storage$subject_id)

################################################################################

# edit flag to note mismatches/instances where leading zeros were re-introduced (CBR)
manifest_storage$flag <- ifelse(is.na(manifest_storage$flag), "", manifest_storage$flag)
manifest_storage$flag <- ifelse(manifest_storage$received_source == "CBR" & manifest_storage$subject_id_length < 9, 
                                paste0(manifest_storage$flag, " ", "MRN < 9 digits + leading 0s restored"), manifest_storage$flag)
manifest_storage$flag <- trimws(manifest_storage$flag)
manifest_storage$flag <- ifelse(manifest_storage$flag == "", NA, manifest_storage$flag)

# add in those leading zeros in cases (CBR)

manifest_storage$subject_id <- ifelse(manifest_storage$received_source == "CBR" & manifest_storage$subject_id_length < 9, 
                                      with_options(c(scipen = 999), str_pad(manifest_storage$subject_id, 9, pad = "0")), 
                                      manifest_storage$subject_id)

################################################################################

# edit flag to note mismatches/instances where leading zeros were re-introduced (CSTP)
manifest_storage$flag <- ifelse(is.na(manifest_storage$flag), "", manifest_storage$flag)
manifest_storage$flag <- ifelse(manifest_storage$received_source == "CSTP" & manifest_storage$subject_id_length < 8, 
                                paste0(manifest_storage$flag, " ", "UMID < 8 digits + leading 0s restored"), manifest_storage$flag)
manifest_storage$flag <- trimws(manifest_storage$flag)
manifest_storage$flag <- ifelse(manifest_storage$flag == "", NA, manifest_storage$flag)

# add in those leading zeros in cases (CBR)

manifest_storage$subject_id <- ifelse(manifest_storage$received_source == "CSTP" & manifest_storage$subject_id_length < 8, 
                                      with_options(c(scipen = 999), str_pad(manifest_storage$subject_id, 8, pad = "0")), 
                                      manifest_storage$subject_id)

################################################################################

# edit flag to note mismatches/instances where leading zeros were re-introduced (ED_IDNOW)
manifest_storage$flag <- ifelse(is.na(manifest_storage$flag), "", manifest_storage$flag)
manifest_storage$flag <- ifelse(manifest_storage$received_source == "ED_IDNOW" & manifest_storage$subject_id_length < 9, 
                                paste0(manifest_storage$flag, " ", "MRN < 9 digits + leading 0s restored"), manifest_storage$flag)
manifest_storage$flag <- trimws(manifest_storage$flag)
manifest_storage$flag <- ifelse(manifest_storage$flag == "", NA, manifest_storage$flag)

# add in those leading zeros in cases (CBR)

manifest_storage$subject_id <- ifelse(manifest_storage$received_source == "ED_IDNOW" & manifest_storage$subject_id_length < 9, 
                                      with_options(c(scipen = 999), str_pad(manifest_storage$subject_id, 9, pad = "0")), 
                                      manifest_storage$subject_id)


################################################################################
#                           File Write-Outs                                    #
################################################################################

### write compiled manifest file out
### in this case, we'll always overwrite the old file, if it does exist
write.csv(manifest_storage, paste0(outputLOC, "/sample_full_manifest_list.csv"), row.names = FALSE, na = "")

### write output report

### create date formatting

# # add leading zero to month
# if (length(month(Sys.Date()))){
#   m <- paste0("0", month(Sys.Date()))
# } else {
#   m <- month(Sys.Date())
# }
# # add leading zero to day
# if (length(day(Sys.Date()))){
#   d <- paste0("0", day(Sys.Date()))
# } else {
#   d <- day(Sys.Date())
# }

today <- current_date_string()

wb <- loadWorkbook(paste0(outputLOC, "/manifest_output_report_template.xlsx"))

writeData(wb, today, sheet = "SUMMARY", startRow = 1, startCol = 2)

#writeData(wb, each_file_row_count, sheet = "SUMMARY", startRow = 3, startCol = 2)
writeData(wb, nrow(manifest_storage), sheet = "SUMMARY", startRow = 4, startCol = 2)

writeData(wb, ncol(manifest_storage), sheet = "SUMMARY", startRow = 6, startCol = 2)

if (exists("duplicate_ssc")){
  writeData(wb, duplicate_ssc, sheet = "DUPLICATES",   
                 startRow = 1, startCol = 1)
}

zeros <- rbind(filter(manifest_storage, grepl("MRN < 9", flag)), filter(manifest_storage, grepl("UMID < 8", flag)))
writeData(wb, zeros, sheet = "RESTORE_ZEROS",   
          startRow = 1, startCol = 1)

saveWorkbook(wb, paste0(outputLOC, "/manifest_output_report_", today, ".xlsx"), overwrite = TRUE)
