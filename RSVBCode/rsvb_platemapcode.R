################################################################################
#       Creation of Plate Map Dataset for Influenza Genetic Sampling           #
#                         Created: 11/17/2021                                  #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
library(openxlsx)

options(scipen=999)

################################################################################
#                Plate Map Files - Upload and Data Checks                      #
################################################################################

# plate map file path
platemap_fp <- paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/PlateMaps")

### output location of plate map files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/PlateMaps/PlateMapsComplete")

################################################################################


file_list <- list.files(pattern = "*.xlsx", path = platemap_fp)

plate_storage <- data.frame()

for (each_page in file_list){
  plate1 <- read.xlsx(paste0(platemap_fp, "/", each_page), sheet = 1)
  
  plate1 <- plate1[, 1:7]
  plate_storage <- rbind(plate_storage, plate1)
}

#### Split up processing plate column

### get plate creation date
plate_storage$Plate_Date <- paste0(substr(plate_storage$Processing.Plate, 1, 4), "-", substr(plate_storage$Processing.Plate, 5, 6), "-", substr(plate_storage$Processing.Plate, 7, 8))
### get plate run platform
plate_storage$Plate_Platform <- sapply(strsplit(as.character(plate_storage$Processing.Plate),'_'), "[", 3)
### get plate number
plate_storage$Plate_Number <- sapply(strsplit(as.character(plate_storage$Processing.Plate),'_'), "[", 5)

### Source Formatting
### Source should refer to the manifest report that the sample came in with
# remove commas, if they are there
plate_storage$Source <- gsub(",", " ", plate_storage$Source)

### check Source column for following format: <character string>, <space>, date as M-D-YYYY
space_check <- grepl(" ", plate_storage$Source)
if (any(space_check == FALSE)){
  print("Warning: There are some Plate records without space character separater in Source column.")
}

date_part <- sapply(strsplit(as.character(plate_storage$Source), " "), "[", 2)
date_part <- as.POSIXct(date_part, format = "%m-%d-%Y")
if (any(is.na(date_part))){
  print("Warning: There are some Plate records without date information in Source column.")
}

#### Get source information to use for joining with the manifest record
plate_storage$Source_Date <- sapply(strsplit(as.character(plate_storage$Source),' '), "[", 2)
plate_storage$Source_Date <- as.character(as_date(paste0(sapply(strsplit(as.character(plate_storage$Source_Date),'-'), "[", 3), "-", sapply(strsplit(as.character(plate_storage$Source_Date),'-'), "[", 1), "-", sapply(strsplit(as.character(plate_storage$Source_Date),'-'), "[", 2))))

plate_storage$Source_Location <- sapply(strsplit(as.character(plate_storage$Source),' '), "[", 1)

plate_storage$Source_Location <- ifelse(is.na(plate_storage$Source_Date), plate_storage$Source, plate_storage$Source_Location)

# note, do not need to add leading zeros, subject ID is not used here.

### drop columns
plate_storage <- plate_storage %>% select(Processing.Plate, Slot, Sample.ACCN, Barcode, Source_Date, Source_Location, Plate_Date, Plate_Platform, Plate_Number)

### rename columns 
colnames(plate_storage) <- c("PlateName", "PlatePosition", "SampleID", "SampleBarcode", "SampleSourceDate", "SampleSourceLocation", "PlateDate", "PlatePlatform", "PlateNumber")

### write out compiled plate information
write.csv(plate_storage, paste0(outputLOC, "/sample_full_plate_list.csv"), row.names = FALSE, na = "")
