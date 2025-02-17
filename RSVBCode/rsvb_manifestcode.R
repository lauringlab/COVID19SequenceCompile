################################################################################
#       Creation of Manifest Dataset for Influenza A Genetic Sampling          #
#                          Created: 11/17/2021                                 #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
library(openxlsx)
library(withr)

################################################################################
#                Manifest Files - Upload and Data Checks                       #
################################################################################

# Manifest file paths (there should be a path per source)
ivy4_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/CDCIVY4")
ivy5_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/CDCIVY5")
ivy6_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/CDCIVY6")
ivy7_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/CDCIVY7")
hive_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/HIVE")
right_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/RIGHT")
martin_manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/MARTIN")

manifest_folder_list <- c(ivy4_manifest_fp, ivy5_manifest_fp, ivy6_manifest_fp, ivy7_manifest_fp, hive_manifest_fp, right_manifest_fp, martin_manifest_fp)

### output location of manifest files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/ManifestsComplete")

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
  
  if (each_folder == ivy4_manifest_fp){
    # process ivy manifest
    
    # read in manifests
    file_list <- list.files(pattern = "*.xlsx", path = each_folder)
    for (ivym in file_list){
      file_in <- read.xlsx(paste0(each_folder, "/", ivym), detectDates = TRUE)
      
      # "position", "sample_id", "subject_id", "coll_date", "flag"
      file_in <- file_in %>% select(`Position.#`, Aliquot.ID, `Study.ID`, `Collection.Date`, Comments)
      colnames(file_in) <- c("position", "sample_id", "subject_id", "coll_date", "flag")
      
      # sometimes have to cut time off of collection date (from excel)
      file_in$coll_date <- substr(as.character(file_in$coll_date), 1, 10)
      
      # add in 2 new columns: received_date and received_source (from file name)
      file_in$received_date <- date_from_file(ivym)
      
      
      rec_source <- trimws(as.character(strsplit(ivym, "_")[[1]][1]))
      file_in$received_source <- rec_source
      file_in$coll_date <- as.character(file_in$coll_date)
      
      # bind all rows together
      manifest_storage <- rbind(manifest_storage, file_in)
      
    }
  } else if (each_folder == ivy5_manifest_fp){
        # process ivy manifest
        
        # read in manifests
        file_list <- list.files(pattern = "*.xlsx", path = each_folder)
        for (ivym in file_list){
          file_in <- read.xlsx(paste0(each_folder, "/", ivym), detectDates = TRUE)
          
          # "position", "sample_id", "subject_id", "coll_date", "flag"
          file_in <- file_in %>% select(`Position.#`, Aliquot.ID, `Study.ID`, `Collection.Date`, Comments)
          colnames(file_in) <- c("position", "sample_id", "subject_id", "coll_date", "flag")
          
          # sometimes have to cut time off of collection date (from excel)
          file_in$coll_date <- substr(as.character(file_in$coll_date), 1, 10)
          
          # add in 2 new columns: received_date and received_source (from file name)
          file_in$received_date <- date_from_file(ivym)
          
          
          rec_source <- trimws(as.character(strsplit(ivym, "_")[[1]][1]))
          file_in$received_source <- rec_source
          file_in$coll_date <- as.character(file_in$coll_date)
          
          # bind all rows together
          manifest_storage <- rbind(manifest_storage, file_in)
          
        }
  } else if (each_folder == ivy6_manifest_fp){
    # process ivy manifest
    
    # read in manifests
    file_list <- list.files(pattern = "*.xlsx", path = each_folder)
    for (ivym in file_list){
      file_in <- read.xlsx(paste0(each_folder, "/", ivym), detectDates = TRUE)
      
      # "position", "sample_id", "subject_id", "coll_date", "flag"
      file_in <- file_in %>% select(`Position.#`, Aliquot.ID, `Study.ID`, `Collection.Date`, Comments)
      colnames(file_in) <- c("position", "sample_id", "subject_id", "coll_date", "flag")
      
      # sometimes have to cut time off of collection date (from excel)
      file_in$coll_date <- substr(as.character(file_in$coll_date), 1, 10)
      
      # add in 2 new columns: received_date and received_source (from file name)
      file_in$received_date <- date_from_file(ivym)
      
      
      rec_source <- trimws(as.character(strsplit(ivym, "_")[[1]][1]))
      file_in$received_source <- rec_source
      file_in$coll_date <- as.character(file_in$coll_date)
      
      # bind all rows together
      manifest_storage <- rbind(manifest_storage, file_in)
      
    }

  } else if (each_folder == ivy7_manifest_fp){
    # process ivy manifest
    
    # read in manifests
    file_list <- list.files(pattern = "*.xlsx", path = each_folder)
    for (ivym in file_list){
      file_in <- read.xlsx(paste0(each_folder, "/", ivym), detectDates = TRUE)
      
      # "position", "sample_id", "subject_id", "coll_date", "flag"
      file_in <- file_in %>% select(`Position.#`, Aliquot.ID, `Study_ID`, `Collection.Date`, Comments)
      colnames(file_in) <- c("position", "sample_id", "subject_id", "coll_date", "flag")
      
      # sometimes have to cut time off of collection date (from excel)
      file_in$coll_date <- substr(as.character(file_in$coll_date), 1, 10)
      
      # add in 2 new columns: received_date and received_source (from file name)
      file_in$received_date <- date_from_file(ivym)
      
      
      rec_source <- trimws(as.character(strsplit(ivym, "_")[[1]][1]))
      file_in$received_source <- rec_source
      file_in$coll_date <- as.character(file_in$coll_date)
      
      # bind all rows together
      manifest_storage <- rbind(manifest_storage, file_in)
      
    }
        
   } else if (each_folder == right_manifest_fp){
      #process RIGHT manifest
     full_right <- data.frame()
     # read in manifests
     file_list <- list.files(pattern = "*.csv", path = each_folder)
     for (rightb in file_list){
       right_one <- read.csv(paste0(each_folder, "/", rightb), colClasses = "character")
       
       rfile_one <- rbind(right_one, full_right)
       
       rfile_one <- rfile_one %>% select(specimen_id, study_id, date_of_collection, site_name,
                                         specimen_type, manifest_creation_date, record_id, rsv_subtype)
       
       # sometimes have to cut time off of collection date (from excel)
       rfile_one$date_of_collection <- substr(as.character(rfile_one$date_of_collection), 1, 10)
       
       
       colnames(rfile_one) <- c("sample_id", "subject_id", "coll_date", "site",
                                "specimen_type", "received_date", "record_id", "flag")
       # add in 2 new columns: position and received_source (from file name)
       rfile_one$position <- ""
       
       
       rec_source <- trimws(as.character(strsplit(rightb, "_")[[1]][1]))
       rfile_one$received_source <- rec_source
       rfile_one$coll_date <- as.character(rfile_one$coll_date)
       
       
       #select the columns in the correct order
       rfile_one <- rfile_one %>% select(position, sample_id, subject_id, coll_date, flag, received_date, received_source)
       
       # bind all rows together
       manifest_storage <- rbind(manifest_storage, rfile_one)
      }
    
  } else {
  
  ### get names of all .csv files in folder
  file_list <- list.files(pattern = "*.csv", path = each_folder)
  
    if (length(file_list) != 0){
      
          # then iterate through files within each folder
          for (each_file in file_list){
            #print(each_file)
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
            
          
            ## reformat coll_date to YYYY-MM-DD format if necessary
            test_date_format <- substr(as.character(file_in[1, 4]), 1, 4)
            #print(test_date_format)
            
            if (is.na(as.numeric(test_date_format))){
              #a <- file_in$coll_date
              file_in$coll_date <- as.POSIXct(file_in$coll_date, format = "%m/%d/%y")
            }
            
            if(any(is.na(file_in$coll_date))){
              #print(each_file)
              #print("There are missing collection dates.")
              
              ## fill the missings with current date, and edit flag
              file_in <- file_in %>% mutate(flag = case_when(is.na(coll_date) ~ gsub("NA", "", paste0(flag, "Missing Date in Manifest - Replaced with Today Date")), 
                                                             T ~ flag), 
                                            coll_date = case_when(is.na(coll_date) ~ as.character(Sys.Date()), 
                                                                  T ~ as.character(coll_date)))
            }
            
            # add in 2 new columns: received_date and received_source (from file name)
            file_in$received_date <- date_from_file(each_file)
            
            rec_source <- trimws(as.character(strsplit(each_file, "_")[[1]][1]))
            file_in$received_source <- rec_source
            
            file_in$coll_date <- as.character(file_in$coll_date)
            
            # bind all rows together
            manifest_storage <- rbind(manifest_storage, file_in)
          
          }
          
        
          ### select only distinct rows
          manifest_storage <- manifest_storage %>% distinct()
    } else {
      print(paste0("No files in folder = ", each_folder))
    }
  }
}

manifest_storage <- filter(manifest_storage, !is.na(sample_id) & !is.na(subject_id))
manifest_storage$coll_date <- as.character(manifest_storage$coll_date)


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


####################
# collection date formatting

# manifest_storage <- manifest_storage %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
#                                           grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
#                                           T ~ NA_character_))

################################################################################
# check for instances where leading zeros may have been dropped from subject id
# CBR should == 9 digits (MRN)
# ED NOW should == 9 digits (MRN)
# CSTP should == 8 digits (UMID)

# add a column to check length of subject_id
manifest_storage$subject_id_length <- nchar(manifest_storage$subject_id)

manifest_storage <- subject_id_length_QA(manifest_storage, "CBR")

################################################################################
#                           File Write-Outs                                    #
################################################################################

### write compiled manifest file out
### in this case, we'll always overwrite the old file, if it does exist
write.csv(manifest_storage, paste0(outputLOC, "/sample_full_manifest_list.csv"), row.names = FALSE, na = "")

### write output report

### create date formatting
# today <- current_date_string()
# 
# wb <- loadWorkbook(paste0(outputLOC, "/manifest_output_report_template.xlsx"))
# 
# writeData(wb, today, sheet = "SUMMARY", startRow = 1, startCol = 2)
# 
# #writeData(wb, each_file_row_count, sheet = "SUMMARY", startRow = 3, startCol = 2)
# writeData(wb, nrow(manifest_storage), sheet = "SUMMARY", startRow = 4, startCol = 2)
# 
# writeData(wb, ncol(manifest_storage), sheet = "SUMMARY", startRow = 6, startCol = 2)
# 
# if (exists("duplicate_ssc")){
#   writeData(wb, duplicate_ssc, sheet = "DUPLICATES",   
#             startRow = 1, startCol = 1)
# }
# 
# zeros <- rbind(filter(manifest_storage, grepl("MRN < 9", flag)), filter(manifest_storage, grepl("UMID < 8", flag)))
# writeData(wb, zeros, sheet = "RESTORE_ZEROS",   
#           startRow = 1, startCol = 1)
# 
# miss_dats <- filter(manifest_storage, grepl("Missing Date", flag))
# writeData(wb, miss_dats, sheet = "MISSING_DATES", startRow = 1, startCol = 1)
# 
# saveWorkbook(wb, paste0(outputLOC, "/manifest_output_report_", today, ".xlsx"), overwrite = TRUE)
