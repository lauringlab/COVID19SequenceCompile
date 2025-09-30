################################################################################
#         Creation of Manifest Dataset for hMPV Genetic Sampling           #
#                           Created: 11/4/2024                                #
#                 Code Edited By: Leigh Papalambros                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
library(openxlsx)
library(withr)
library(stringr)


################################################################################
#                Manifest Files - Upload and Data Checks                       #
################################################################################

#starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"

# Manifest file paths (there should be a path per source)
cdcivy_manifest_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/Manifests/CDCIVY")
cdcivy7_manifest_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/Manifests/CDCIVY7")
hive_manifest_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/Manifests/HIVE")

### output location of manifest files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/Manifests/ManifestsComplete")

# each manifest provider will have their own folder, with all files inside

################################################################################

# read in list of manifest files that have already been processed in prev file
#processed_manifest_file_names <- readRDS(paste0(outputLOC, "/current_manifest_list.RDS"))
#processed_manifest_file_names <- Filter(function(x) !any(grepl("TRIN", x)), processed_manifest_file_names)

################################################################################


manifest_storage <- data.frame()

################################################################################
## handle cdc ivy manifests

# read in cdc ivy site code list 
cdc_sites <- read.csv(paste0(cdcivy_manifest_fp, "/Keys/CDC_SiteCodebook.csv"), colClasses = "character")

#cdc_file_list <- list.files(pattern = "*.xlsx", path = cdcivy_manifest_fp)

cdc_file_list24 <- list.files(pattern = "*.xlsx", path = cdcivy_manifest_fp)
cdc_file_list25 <- list.files(pattern = "*.xlsx", path = cdcivy7_manifest_fp)
#print(cdc_file_list24)

print("Processing CDCIVY Manifests")

cdc_file_list_24 <- c()
cdc_file_list_25 <- c()

#print(cdc_file_list24)
#for (each_file in cdc_file_list24){
#  if (each_file %in% processed_manifest_file_names){
#    xx <- "skip it"
#  } else {
#    cdc_file_list <- c(cdc_file_list, each_file)
#  }
#}


cdc_ivy_storage <- data.frame()
fileone <- data.frame()
filetwo <- data.frame()
full_ivy <- data.frame()
full_ivy_2 <- data.frame()
all_ivy <- data.frame()


################################################################################

## processing ivy5 and ivy6 manifests ##
for (each_file in cdc_file_list24){
  fileone <- read.xlsx(paste0(cdcivy_manifest_fp, "/", each_file), sheet = 1, detectDates = TRUE)
  rec_date <- trimws(as.character(strsplit(each_file, "_")[[1]][3]))
  fileone$received_date <- rec_date
  rec_source <- trimws(as.character(strsplit(each_file, "_")[[1]][1]))
  fileone$received_source <- rec_source
  fileone <- filter(fileone, !is.na(as.numeric(`Position.#`)))
  
  if (any(!grepl("-", fileone$Collection.Date))){
    fileone <- fileone %>% mutate(Collection.Date = case_when(grepl("/", Collection.Date) ~ as.character(as.POSIXct(Collection.Date, format = "%m/%d/%Y")), 
                                                              !grepl("/", Collection.Date) & !grepl("-", Collection.Date) ~ as.character(as.Date(Collection.Date, origin = "1899-12-30")), 
                                                              T ~ Collection.Date))
  }
  
  fileone$Collection.Date <- as.character(paste0(substr(fileone$Collection.Date, 1, 4), "-", substr(fileone$Collection.Date, 6, 7), "-", substr(fileone$Collection.Date, 9, 10)))
  
  ## keep full set of cdc ivy rows separate, for back checks on full data
  full_ivy <- rbind(full_ivy, fileone)
  ## change all column names to lowercase and remove leading/lagging white space
  ## to make it easier to process
  cdc_names <- colnames(fileone)
  new_names <- c()
  for (i in cdc_names){
    new_names <- c(new_names, tolower(trimws(i)))
  }
  colnames(fileone) <- new_names
  
  fileone <- full_ivy %>% select(`Position.#`, Study.ID, Collection.Date, Aliquot.ID, received_source, received_date)
  
  colnames(fileone) <- c("position", "subject_id", "coll_date", "sample_id","received_source", "received_date")
  fileone$subject_id <- trimws(fileone$subject_id)
  fileone$subject_id <- as.numeric(fileone$subject_id)
  fileone$coll_date <- as.character(fileone$coll_date)
  
  #fileone <- fileone %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
  #                                                      T ~ as.numeric(substr(subject_id, 3, 4))))
  
  #check_site_codes <- fileone %>% select(subject_id, sample_id) %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
  #                                                                                                 T ~ as.numeric(substr(subject_id, 3, 4))))
  #check_site_codes <- merge(check_site_codes, cdc_sites, by.x = c("site_number"), by.y = c("Number"), all.x = TRUE)
  
 # if(any(is.na(check_site_codes$Institution))){
  #  message(each_file)
 #   message(filter(check_site_codes, is.na(Institution)))
 #   stop("No Site Numerical Match")
  #}
  
  #site_bit <- cdc_sites %>% select(Number, SiteCode)
  #colnames(site_bit) <- c("Number", "SiteName")
  #fileone <- merge(fileone, site_bit, by.x = c("site_number"), by.y = c("Number"))
  
  #fileone$received_date <- date_from_file(cdc_file_list24)
  #fileone$received_date <- str_extract(cdc_file_list24, "\\d{8}")
  
  
  ### add in "regular" manifest columns
  fileone$flag <- NA
  
  ### re-arrange variables
  fileone <- fileone %>% select(position, sample_id, subject_id, coll_date, flag, received_date, received_source)
}

## processing ivy7 manifests ##
for (each_file in cdc_file_list25){
  filetwo <- read.xlsx(paste0(cdcivy7_manifest_fp, "/", each_file), sheet = 1, detectDates = TRUE)
  rec_date <- trimws(as.character(strsplit(each_file, "_")[[1]][3]))
  filetwo$received_date <- rec_date
  rec_source <- trimws(as.character(strsplit(each_file, "_")[[1]][1]))
  filetwo$received_source <- rec_source
  filetwo <- filter(filetwo, !is.na(as.numeric(`Position.#`)))
  
  if (any(!grepl("-", filetwo$Collection.Date))){
    fileone <- filetwo %>% mutate(Collection.Date = case_when(grepl("/", Collection.Date) ~ as.character(as.POSIXct(Collection.Date, format = "%m/%d/%Y")), 
                                                              !grepl("/", Collection.Date) & !grepl("-", Collection.Date) ~ as.character(as.Date(Collection.Date, origin = "1899-12-30")), 
                                                              T ~ Collection.Date))
  }
  
  filetwo$Collection.Date <- as.character(paste0(substr(filetwo$Collection.Date, 1, 4), "-", substr(filetwo$Collection.Date, 6, 7), "-", substr(filetwo$Collection.Date, 9, 10)))
  
  ## keep full set of cdc ivy rows separate, for back checks on full data
  full_ivy_2 <- rbind(full_ivy_2, filetwo)
  ## change all column names to lowercase and remove leading/lagging white space
  ## to make it easier to process
  cdc_names <- colnames(filetwo)
  new_names <- c()
  for (i in cdc_names){
    new_names <- c(new_names, tolower(trimws(i)))
  }
  colnames(filetwo) <- new_names
  
  filetwo <- full_ivy_2 %>% select(`Position.#`, Study_ID, Collection.Date, Aliquot.ID, received_source, received_date)
  
  colnames(filetwo) <- c("position", "subject_id", "coll_date", "sample_id","received_source", "received_date")
  filetwo$subject_id <- trimws(filetwo$subject_id)
  filetwo$subject_id <- as.numeric(filetwo$subject_id)
  filetwo$coll_date <- as.character(filetwo$coll_date)
  
  #filetwo <- filetwo %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
  #                                                      T ~ as.numeric(substr(subject_id, 3, 4))))
  
  #check_site_codes <- filetwo %>% select(subject_id, sample_id) %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
  #                                                                                                 T ~ as.numeric(substr(subject_id, 3, 4))))
  #check_site_codes <- merge(check_site_codes, cdc_sites, by.x = c("site_number"), by.y = c("Number"), all.x = TRUE)
  
  #if(any(is.na(check_site_codes$Institution))){
  #  message(each_file)
  #  message(filter(check_site_codes, is.na(Institution)))
 #   stop("No Site Numerical Match")
 # }
  
 # site_bit <- cdc_sites %>% select(Number, SiteCode)
 # colnames(site_bit) <- c("Number", "SiteName")
 # filetwo <- merge(filetwo, site_bit, by.x = c("site_number"), by.y = c("Number"))
  
  #fileone$received_date <- date_from_file(cdc_file_list24)
  #fileone$received_date <- str_extract(cdc_file_list24, "\\d{8}")
  
  
  ### add in "regular" manifest columns
  filetwo$flag <- NA
  
  ### re-arrange variables
  filetwo <- filetwo %>% select(position, sample_id, subject_id, coll_date, flag, received_date, received_source)
}


## selecting the files that are needed from the fully IVY files

fullivy_filt <- full_ivy %>% select('Position.#', Study.ID, Aliquot.ID, Collection.Date, Comments, received_source, received_date)
colnames(fullivy_filt) <- c("position", "subject_id", "sample_id", "coll_date", "flag", "received_source", "received_date")

fullivy2_filt <- full_ivy_2 %>% select('Position.#', Study_ID, Aliquot.ID, Collection.Date, Comments, received_source, received_date)
colnames(fullivy2_filt) <- c("position", "subject_id", "sample_id", "coll_date","flag", "received_source", "received_date")

all_ivy_wf <- rbind(fullivy_filt, fullivy2_filt)

## format the received_date into date format
#str(all_ivy_wf)
all_ivy_wf$received_date <- as.Date(all_ivy_wf$received_date, format = "%Y%m%d")

all_ivy_wf <- all_ivy_wf %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
                                                        T ~ as.numeric(substr(subject_id, 3, 4))))
  
check_site_codes <- all_ivy_wf %>% select(subject_id, sample_id) %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
                                                                                                   T ~ as.numeric(substr(subject_id, 3, 4))))
check_site_codes <- merge(check_site_codes, cdc_sites, by.x = c("site_number"), by.y = c("Number"), all.x = TRUE)
  
if(any(is.na(check_site_codes$Institution))){
    message(each_file)
    message(filter(check_site_codes, is.na(Institution)))
    stop("No Site Numerical Match")
  }
  
  site_bit <- cdc_sites %>% select(Number, SiteCode)
  colnames(site_bit) <- c("Number", "SiteName")
  all_ivy_wf <- merge(all_ivy_wf, site_bit, by.x = c("site_number"), by.y = c("Number"))

wf_ivy <- all_ivy_wf %>% select("position", "subject_id", "sample_id", "coll_date", "flag", "received_source", "received_date", "SiteName")
  
  
  ### HIVE manifest section ###

hive_file_list <- c()

hive_file_list <- list.files(pattern = "*.xlsx", path = hive_manifest_fp)
#print(cdc_file_list24)

print("Processing HIVE Manifests")



full_hive <- data.frame()
hivefile <- data.frame()


for (each_hive in hive_file_list){
  hivefile <- read.xlsx(paste0(hive_manifest_fp, "/", each_hive), sheet = 1, detectDates = TRUE)
  
  hivefile$SPEC_COLLECT_DATE <- as.Date(hivefile$SPEC_COLLECT_DATE, origin = "1899-12-30")
  hivefile$received_date <- ""
  
  #hivefile$received_date <- date_from_file(each_file)
  
  rec_source <- trimws(as.character(strsplit(each_hive, "_")[[1]][1]))
  hivefile$received_source <- rec_source
  hivefile$received_date <- ""
  
  full_hive <- rbind(hivefile, full_hive)
  
  full_hive_wf <- full_hive %>% select(STUDYID, SPEC_COLLECT_DATE, ACCN, SEASON, received_date, received_source)
  full_hive_wf$position <- ""
  full_hive_wf$SiteName <- ""
  full_hive_wf <- full_hive_wf %>% select(position, ACCN, STUDYID, SPEC_COLLECT_DATE, SEASON, received_date, received_source, SiteName)
  
  colnames(full_hive_wf) <- c("position", "sample_id", "subject_id", "coll_date",  "flag", "received_date", "received_source", "SiteName")
  
  
  
}

#str(full_hive_wf)
## formatting the coll_date column to make sure it stays in date format
full_hive_wf$coll_date <- as.character(as.Date(full_hive_wf$coll_date))

### Merge all manifest files into one big file ###

class(full_hive_wf$coll_date)
class(wf_ivy$coll_date)

manifest_storage <- data.frame()

## selecting and formatting IVY columns for merge with other manifest

#m_all_ivy <- all_ivy %>% select('position.#', aliquot.id, collection.date, study.id, comments)
#m_all_ivy$
#colnames(m_all_ivy) <- c("position", "subject_id", "coll_date", "sample_id", "flag", "received_date", "received_source")


manifest_storage <- rbind(wf_ivy, full_hive_wf)


################################################################################
# check for sample_id/subject_id/coll_date duplicates
# count of unique sample_id, subject_id, coll_date combinations

if (nrow(manifest_storage) > 0){
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
}


### write compiled manifest file out
### in this case, we'll always overwrite the old file, if it does exist
write.csv(manifest_storage, paste0(outputLOC, "/sample_full_manifest_list.csv"), row.names = FALSE, na = "")









