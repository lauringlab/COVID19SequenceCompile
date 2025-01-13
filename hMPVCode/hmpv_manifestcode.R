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

################################################################################
#                Manifest Files - Upload and Data Checks                       #
################################################################################

#starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"

# Manifest file paths (there should be a path per source)
cdcivy_manifest_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/Manifests/CDCIVY")


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
#print(cdc_file_list24)

print("Processing CDCIVY Manifests")

cdc_file_list <- c()
#for (each_file in cdc_file_list24){
#  if (each_file %in% processed_manifest_file_names){
#    xx <- "skip it"
#  } else {
#    cdc_file_list <- c(cdc_file_list, each_file)
#  }
#}


cdc_ivy_storage <- data.frame()
full_ivy <- data.frame()


for (each_file in cdc_file_list24){
  #print(each_file)
  if (grepl("CDCIVY5", each_file) | grepl("CDCIVY6", each_file) | grepl("CDCIVY7", each_file)){
    #print("IVY4")
    #print(each_file)
    fileone <- read.xlsx(paste0(cdcivy_manifest_fp, "/", each_file), sheet = 1, detectDates = TRUE)
    
    fileone <- filter(fileone, !is.na(as.numeric(`Position.#`)))
    
    if (any(!grepl("-", fileone$Collection.Date))){
      fileone <- fileone %>% mutate(Collection.Date = case_when(grepl("/", Collection.Date) ~ as.character(as.POSIXct(Collection.Date, format = "%m/%d/%Y")), 
                                                                !grepl("/", Collection.Date) & !grepl("-", Collection.Date) ~ as.character(as.Date(Collection.Date, origin = "1899-12-30")), 
                                                                T ~ Collection.Date))
    }
    
    fileone$Collection.Date <- as.character(paste0(substr(fileone$Collection.Date, 1, 4), "-", substr(fileone$Collection.Date, 6, 7), "-", substr(fileone$Collection.Date, 9, 10)))
    
    
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
    
    fileone <- fileone %>% select(`position.#`, study.id, collection.date, aliquot.id)
    
    colnames(fileone) <- c("position", "subject_id", "coll_date", "sample_id")
    fileone$subject_id <- trimws(fileone$subject_id)
    fileone$subject_id <- as.numeric(fileone$subject_id)
    fileone$coll_date <- as.character(fileone$coll_date)
    
    fileone <- fileone %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
                                                          T ~ as.numeric(substr(subject_id, 3, 4))))
    
    check_site_codes <- fileone %>% select(subject_id, sample_id) %>% mutate(site_number = case_when(substr(subject_id, 1, 1) == "C" ~ as.numeric(substr(subject_id, 4, 5)), 
                                                                                                     T ~ as.numeric(substr(subject_id, 3, 4))))
    check_site_codes <- merge(check_site_codes, cdc_sites, by.x = c("site_number"), by.y = c("Number"), all.x = TRUE)
    
    if(any(is.na(check_site_codes$Institution))){
      message(each_file)
      message(filter(check_site_codes, is.na(Institution)))
      stop("No Site Numerical Match")
    }
    
    site_bit <- cdc_sites %>% select(Number, SiteCode)
    colnames(site_bit) <- c("Number", "SiteName")
    fileone <- merge(fileone, site_bit, by.x = c("site_number"), by.y = c("Number"))
    
    fileone$received_date <- date_from_file(each_file)
    
    rec_source <- trimws(as.character(strsplit(each_file, "_")[[1]][1]))
    fileone$received_source <- rec_source
    
    ### add in "regular" manifest columns
    fileone$flag <- NA
    
    ### re-arrange variables
    fileone <- fileone %>% select(position, sample_id, subject_id, coll_date, flag, received_date, received_source, SiteName)
    
  } else {
  
  cdc_ivy_storage <- rbind(cdc_ivy_storage, fileone)
  }
}

### write out full ivy set
write.csv(full_ivy, paste0(cdcivy_manifest_fp, "/Full_IVY_Set/IVY_sample_full_manifest_list24.csv"), row.names = FALSE, na = "")
#write.csv(full_ivy4, paste0(cdcivy_manifest_fp, "/Full_IVY_Set/IVY4_sample_full_manifest_list22.csv"), row.names = FALSE, na = "")

### add onto main manifest file HERE
if (nrow(manifest_storage) > 0){
  manifest_storage$SiteName <- NA
}
manifest_storage <- rbind(manifest_storage, cdc_ivy_storage)

### write compiled manifest file out
### in this case, we'll always overwrite the old file, if it does exist
write.csv(manifest_storage, paste0(outputLOC, "/sample_full_manifest_list.csv"), row.names = FALSE, na = "")









