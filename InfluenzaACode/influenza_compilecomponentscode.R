################################################################################
#       Creation of Cumulative Dataset for Influenza Genetic Sampling          #
#                            Created: 2021-11-19                               #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(withr)

################################################################################
#                 Component Files - Upload and Data Checks                     #
################################################################################

# manifest file path
manifest_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/ManifestsComplete")
# platemap file path
plate_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/PlateMaps/PlateMapsComplete")
# nextclade file path
nc_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
# gisaid file path
#gisaid_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
# previous 2021 file path
#prev_2021 <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PreviousLists")

### output location for files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary")

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
nextclade <- read.csv(paste0(nc_fp, "/sample_full_nextclade_list.csv"), colClasses = "character")
#gisaid <- read.csv(paste0(gisaid_fp, "/sample_full_gisaid_list.csv"), colClasses = "character")

#mani_plate_pang_g <- merge(mani_plate_pang, gisaid, by.x = c("sample_id"), by.y = c("sample_id"), all.x = TRUE)
mppnc <- merge(mani_plate, nextclade, by.x = c("sample_id"), by.y = c("SampleID"), all.x = TRUE)

### add column for time in days from plate to nextclade
mppnc$PlateToNextclade_days <- difftime(mppnc$nextclade_HA_runDate, mppnc$PlateDate, units = "days")

################################################################################
## Additional subject_id length check (leading zeros)
## CSTP == 8 (UMIDs), CBR == 9 (MRNs)

# add a column to check length of subject_id
mppnc$subject_id_length <- nchar(mppnc$subject_id)

mppnc <- subject_id_length_QA(mppnc, "CBR")

################################################################################
## add a column to number multiple sample_ids per subject_id

mppnc <- mppnc %>% group_by(subject_id) %>% arrange(coll_date) %>% mutate(sample_per_subject = row_number())

################################################################################
## apply logic for mis-matched pangolin/nextclade info
## necessary for instances where a sample was run on two different plates
## unsure how this logic will work for flu, since not everything gets a nextclade entry

colnames(mppnc)


mppnc2 <- mppnc %>% select(sample_id, subject_id, coll_date,                   
                           flag, received_source, received_date, SampleBarcode,               
                           PlateDate, PlatePlatform, PlateNumber,
                           nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,      
                           nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                           nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                           nextclade_HA_runDate, nextclade_HA_type, 
                           #gisaid_strain, gisaid_epi_isl,              
                           subject_id_length, position, PlateName, PlatePosition,               
                           SampleSourceLocation, PlateToNextclade_days, 
                           sample_per_subject)

################################################################################
### negative control well warning
### this doesn't work with nextclade system

# neg_control <- unique(filter(mppnc2, grepl("NC_", sample_id) & is.na(as.numeric(sample_id)))$sample_id)
# helas <- unique(filter(mppnc2, grepl("hela", tolower(sample_id)))$sample_id)
# 
# check_NCs <- filter(mppnc2, sample_id %in% neg_control | sample_id %in% helas)
# 
# # We want to make sure with each plate that the three negative controls have ???10% of genome covered. 
# 
# check_NCs$neg_control_warning <- ifelse(check_NCs$nextclade_completeness >= 10, 1, 0)
# 
# keep_NCs <- table(check_NCs$PlateName, check_NCs$neg_control_warning)
# 
# write.table(keep_NCs, paste0(outputLOC, "/ReportNotifications/negative_control_warnings.tsv"), sep = "\t")

################################################################################

write.csv(mppnc2, paste0(outputLOC, "/full_compiled_data.csv"), row.names = FALSE, na = "")
write.csv(mppnc2, paste0(outputLOC, "/secret/full_compiled_data.csv"), row.names = FALSE, na = "")
