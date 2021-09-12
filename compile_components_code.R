################################################################################
#       Creation of Cumulative Dataset for COVID-19 Genetic Sampling           #
#                         Last Updated: 5 August 2021                          #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(withr)

################################################################################
#                 Component Files - Upload and Data Checks                     #
################################################################################

# manifest file path
manifest_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/ManifestsComplete")
# platemap file path
plate_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/PlateMapsComplete")
# nextclade file path
nc_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
# pangolin file path
pang_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
# gisaid file path
gisaid_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
# previous 2021 file path
prev_2021 <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PreviousLists")

### output location for files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary")

################################################################################


# first, read in manifest file
manifest <- read.csv(paste0(manifest_fp, "/sample_full_manifest_list.csv"), colClasses = "character")

# read in plate map file
plate_map <- read.csv(paste0(plate_fp, "/sample_full_plate_list.csv"), colClasses = "character")
plate_map <- filter(plate_map, SampleID != "" & !is.na(SampleID)) # remove any rows where sample ID is missing on the plate map
#plate_map$SampleID <- as.character(trimws(plate_map$SampleID))

# store unique number of sample ids 
plate_map_ids <- nrow(plate_map %>% group_by(SampleID, SampleSourceDate) %>% summarize(count = length(SampleID)))


################################################################################
# Warning for if plate map date and manifest date are DIFFERENT

manifest_options <- filter(manifest, received_date != "" & !is.na(received_date)) %>% select(sample_id, subject_id, received_date)
platemap_options <- filter(plate_map, SampleSourceDate != "" & !is.na(SampleSourceDate)) %>% select(SampleID, SampleSourceDate, PlateNumber)

compare_options <- merge(manifest_options, platemap_options, by.x = c("sample_id"), by.y = c("SampleID"))

compare_options$different <- ifelse(compare_options$received_date != compare_options$SampleSourceDate, 1, 0)

compare_options <- filter(compare_options, as_date(received_date) >= as_date("2021-07-01") & different == 1)

if (nrow(compare_options) > 0){
  print(compare_options)
  stop("Mismatched received dates between manifest and platemap.")
}

################################################################################


# merge on plate map file, and only keep rows where plate map file has a sample
mani_plate <- merge(manifest, plate_map, by.x = c("sample_id", "received_date"), by.y = c("SampleID", "SampleSourceDate"), all.y = TRUE)
mani_plate <- filter(mani_plate, !is.na(received_date) & received_date != "")

### double catch in case there is no received date indicated by the Plate Map file
dc <- filter(plate_map, SampleSourceDate == "" | is.na(SampleSourceDate))
mani_plate2 <- merge(manifest, dc, by.x = c("sample_id"), by.y = c("SampleID"), all.y = TRUE)
mani_plate2 <- mani_plate2[ , !names(mani_plate2) %in% c("SampleSourceDate")]

if (nrow(mani_plate2) != nrow(dc)){
  stop("There are duplicate sample_ids between dc and manifest.")
}

mani_plate <- rbind(mani_plate, mani_plate2)

#if (nrow(mani_plate) != plate_map_ids){
#  stop("There are more or less rows in our manifest + plate combination than there were date/sample id combinations in the original plate map file")
#}

#missings <- filter(mani_plate, is.na(subject_id))
#write.csv(missings, "C:/Users/juliegil/Documents/UofM_Work/Lauring_Lab/check_miss_subjects.csv", na = "", row.names = FALSE)

# then, read in pangolin, gisaid, and next clade files
pangolin <- read.csv(paste0(pang_fp, "/sample_full_pangolin_list.csv"), colClasses = "character")
nextclade <- read.csv(paste0(nc_fp, "/sample_full_nextclade_list.csv"), colClasses = "character")
gisaid <- read.csv(paste0(gisaid_fp, "/sample_full_gisaid_list.csv"), colClasses = "character")

# merge these onto mani_plate, always keeping everything from mani_plate
mani_plate_pang <- merge(mani_plate, pangolin, by.x = c("sample_id"), by.y = c("SampleID"), all.x = TRUE)

### add column for time in days from plate to pangolin
mani_plate_pang$PlateToPangolin_days <- difftime(mani_plate_pang$pangolin_runDate, mani_plate_pang$PlateDate, units = "days")

mani_plate_pang_g <- merge(mani_plate_pang, gisaid, by.x = c("sample_id"), by.y = c("sample_id"), all.x = TRUE)
mppnc <- merge(mani_plate_pang_g, nextclade, by.x = c("sample_id"), by.y = c("SampleID"), all.x = TRUE)

### add column for time in days from plate to nextclade
mppnc$PlateToNextclade_days <- difftime(mppnc$nextclade_runDate, mppnc$PlateDate, units = "days")

################################################################################
# create indicator for if Plate to Nextclade or Plate to Pangolin is out of 
# expected range. this will help detect potential "wrong matches" for sample_ids 
# that are sent to us twice & re-run through process

mppnc$IlluminaPangolin_OutOfRange <- ifelse(mppnc$PlatePlatform == "Illumina" & mppnc$PlateToPangolin_days > 8, 1, 0)
mppnc$NanoporePangolin_OutOfRange <- ifelse(mppnc$PlatePlatform == "Nanopore" & mppnc$PlateToPangolin_days > 4, 1, 0)
mppnc$IlluminaNextclade_OutOfRange <- ifelse(mppnc$PlatePlatform == "Illumina" & mppnc$PlateToNextclade_days > 8, 1, 0)
mppnc$NanoporeNextclade_OutOfRange <- ifelse(mppnc$PlatePlatform == "Nanopore" & mppnc$PlateToNextclade_days > 4, 1, 0)

outofrangeset <- filter(mppnc, IlluminaPangolin_OutOfRange == 1 | NanoporePangolin_OutOfRange == 1 | IlluminaNextclade_OutOfRange == 1 | NanoporeNextclade_OutOfRange == 1)

outofrange_output <- filter(mppnc, sample_id %in% outofrangeset$sample_id)

write.csv(outofrange_output, paste0(outputLOC, "/ReportNotifications/out_of_range_alert.csv"), row.names = FALSE, na = "")

################################################################################
## Want to combine with previous 2021 data

prev2 <- read.csv(paste0(prev_2021, "/ProcessedSampleCumulativeList_20210326.csv"), colClasses = "character")

# first, turn MRN & UMID columns into one, sample_id column
prev2$MRN <- gsub("/", "", prev2$MRN)
prev2$MRN <- ifelse(prev2$MRN == "NA", NA, prev2$MRN)
prev2$umid <- gsub("/", "", prev2$umid)
prev2$umid <- ifelse(prev2$umid == "NA", NA, prev2$umid)

prev2$subject_id <- coalesce(prev2$MRN, prev2$umid)

changes <- c("CBR 3-15-2021", "CBR December MFIVE", "Lynx", "LynxDx")
table(prev2$origin)
prev2$changed_origin <- ifelse(prev2$origin %in% changes, 1, 0)
prev2 <- prev2 %>% mutate(note = case_when(changed_origin == 1 ~ paste0(note, "Origin Column changed from ", origin),
                                           T ~ note), 
                          origin = case_when(origin == "CBR 3-15-2021" | origin == "CBR December MFIVE" ~ "CBR", 
                                             origin == "Lynx" | origin == "LynxDx" ~ "CSTP", 
                                             T ~ origin)
                          )


### get plate creation date
prev2$PlateDate <- paste0(substr(prev2$batch, 1, 4), "-", substr(prev2$batch, 5, 6), "-", substr(prev2$batch, 7, 8))

### get plate number
prev2$PlateNumber <- paste0(sapply(strsplit(as.character(prev2$batch),'_'), "[", 3), "_", sapply(strsplit(as.character(prev2$batch),'_'), "[", 4))
prev2$PlateNumber <- ifelse(prev2$PlateNumber == "NA_NA", NA, prev2$PlateNumber)

prev2 <- prev2 %>% select(sample_ID, subject_id, collection_date, note, origin, nanopore_barcode, PlateDate, 
                          platform, 
                          PlateNumber, pangolin_lineage, pangolin_probability, pangolin_status, 
                          pangolin_note, clade, totalMissing, completeness, strain, gisaid_epi_isl)

colnames(prev2) <- c("sample_id", "subject_id", "coll_date", "flag", "received_source", "SampleBarcode", 
                     "PlateDate", "PlatePlatform", "PlateNumber", "pangolin_lineage", "pangolin_probability", 
                     "pangolin_status", "pangolin_note", "nextclade_clade", "nextclade_totalMissing", 
                     "nextclade_completeness", "gisaid_strain", "gisaid_epi_isl")

mppnc2 <- merge(prev2, mppnc, all.x = TRUE, all.y = TRUE)
mppnc2$subject_id <- gsub("/", "", mppnc2$subject_id)
#colnames(mppnc2)

################################################################################
## Additional subject_id length check (leading zeros)
## CSTP == 8 (UMIDs), CBR == 9 (MRNs)

# add a column to check length of subject_id
mppnc2$subject_id_length <- nchar(mppnc2$subject_id)

mppnc2 <- subject_id_length_QA(mppnc2, "CBR")
mppnc2 <- subject_id_length_QA(mppnc2, "EDIDNOW")
mppnc2 <- subject_id_length_QA(mppnc2, "CSTP")

################################################################################
## add a column to number multiple sample_ids per subject_id

mppnc2 <- mppnc2 %>% group_by(subject_id) %>% arrange(coll_date) %>% mutate(sample_per_subject = row_number())

################################################################################

write.csv(mppnc2, paste0(outputLOC, "/full_compiled_data.csv"), row.names = FALSE, na = "")
write.csv(mppnc2, paste0(outputLOC, "/secret/full_compiled_data.csv"), row.names = FALSE, na = "")
