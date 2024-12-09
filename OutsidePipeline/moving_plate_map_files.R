library(tidyverse)
library(openxlsx)

options(scipen=999)

# get plate map from 
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\2_PlateMaps
# and drop it in
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\4_SequenceSampleMetadata\PlateMaps

# remove end code for extended file names
platform <- strsplit(plate_name, "_")[[1]][3]

if (platform == "Illumina"){
  plate_name2 <- paste(strsplit(plate_name, "_")[[1]][1], strsplit(plate_name, "_")[[1]][2], strsplit(plate_name, "_")[[1]][3], strsplit(plate_name, "_")[[1]][4], strsplit(plate_name, "_")[[1]][5], sep = "_")
} else {
  plate_name2 <- plate_name
}

if (grepl("_SC2_", plate_name)){

    file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/2_PlateMaps/", plate_name, ".xlsx"),
                to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), 
              overwrite = TRUE)
  
    # then, read in that excel file that we moved to the inner pipeline set
    file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), sheet = 1)
  
  
} else if (grepl("_IAV_", plate_name)){
  # /Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/2_Plate_Maps
    file.copy(from = paste0(starting_path, "SEQUENCING/INFLUENZA_A/2_Plate_Maps/", plate_name, ".xlsx"),
              to = paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), 
              overwrite = TRUE)
  
    # then, read in that excel file that we moved to the inner pipeline set
    file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), sheet = 1)
  
  
} else if (grepl("_IBV_", plate_name)){
  # /Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/2_Plate_Maps
  file.copy(from = paste0(starting_path, "SEQUENCING/INFLUENZA_B/2_Plate_Maps/", plate_name, ".xlsx"),
            to = paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), 
            overwrite = TRUE)
  
  # then, read in that excel file that we moved to the inner pipeline set
  file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), sheet = 1)
  
  
} else if (grepl("_RSVA_", plate_name)){
  file.copy(from = paste0(starting_path, "SEQUENCING/RSV_A/2_PlateMaps/", plate_name, ".xlsx"),
            to = paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), 
            overwrite = TRUE)
  
  # then, read in that excel file that we moved to the inner pipeline set
  file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), sheet = 1)
  
  
} else if (grepl("_RSVB_", plate_name)){
  file.copy(from = paste0(starting_path, "SEQUENCING/RSV_B/2_PlateMaps/", plate_name, ".xlsx"),
            to = paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), 
            overwrite = TRUE)
  
  # then, read in that excel file that we moved to the inner pipeline set
  file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), sheet = 1)
 
} else if (grepl("_HMPV_", plate_name)){
  file.copy(from = paste0(starting_path, "SEQUENCING/hMPV/2_PlateMaps/", plate_name, ".xlsx"),
            to = paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), 
            overwrite = TRUE)
  
  # then, read in that excel file that we moved to the inner pipeline set
  file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), sheet = 1)
  
} else {
    stop("Plate name does not contain recognized phrase (IAV/SC2).")
}



file_in$Sample.ACCN <- ifelse(grepl('\\.', file_in$Sample.ACCN) & grepl('E', file_in$Sample.ACCN), as.numeric(file_in$Sample.ACCN), file_in$Sample.ACCN)

# remove the well position column
drop_columns <- c("Well.Position")
file_in <- file_in[!names(file_in) %in% drop_columns]

if (grepl("_Illumina_", plate_name)){
  file_in <- file_in[, c(1:6, 10)]
} 

# Rename first 7 columns
colnames(file_in)[1] <- "Processing Plate"
colnames(file_in)[2] <- "Slot"
colnames(file_in)[3] <- "Sample ACCN"
colnames(file_in)[4] <- "Sample MRN"
colnames(file_in)[5] <- "Sample Order#"
colnames(file_in)[6] <- "Barcode"
colnames(file_in)[7] <- "Source"


# need to remove negative control wells from consideration
file_in_source <- filter(file_in, !grepl("control", tolower(file_in$Source)) & !grepl("ctrl", tolower(file_in$Source)))
### check Source column for following format: <character string>, <space>, date as M-D-YYYY
space_check <- grepl(" ", file_in_source$Source)
if (any(space_check == FALSE)){
  print(plate_name)
  print("Warning: There are some records without space character separater in Source column.")
}

date_part <- sapply(strsplit(as.character(file_in_source$Source), " "), "[", 2)
date_part <- as.POSIXct(date_part, format = "%m-%d-%Y")
if (any(is.na(date_part))){
  print(plate_name)
  print("Warning: There are some records without date information in Source column.")
}


if (grepl("_SC2_", plate_name)){
  
  # write the excel file back out
  write.xlsx(file_in, paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), overwrite = TRUE)
  
  
} else if (grepl("_IAV_", plate_name)){
  
  # write the excel file back out
  write.xlsx(file_in, paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), overwrite = TRUE)
  
} else if (grepl("_IBV_", plate_name)){
  
  # write the excel file back out
  write.xlsx(file_in, paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), overwrite = TRUE)
  
} else if (grepl("_RSVA_", plate_name)){
  
  # write the excel file back out
  write.xlsx(file_in, paste0(starting_path, "SEQUENCING/RSV_A/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), overwrite = TRUE)
  
} else if (grepl("_RSVB_", plate_name)){
  
  # write the excel file back out
  write.xlsx(file_in, paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), overwrite = TRUE)
  
} else if (grepl("_HMPV_", plate_name)){
  
  # write the excel file back out
  write.xlsx(file_in, paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/PlateMaps/", plate_name2, ".xlsx"), overwrite = TRUE)
  
} else {
  stop("Plate name does not contain recognized phrase (IAV/SC2/RSV/hMPV).")
}

plate_name <- plate_name2
