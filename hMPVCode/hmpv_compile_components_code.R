################################################################################
#       Creation of Cumulative Dataset for hMPV Genetic Sampling           #
#                         Last Updated: 11/18/2024                         #
#                 Code Edited By: Leigh Papalambros                       #
################################################################################

library(tidyverse)
library(lubridate)
library(withr)
library(arsenal)

################################################################################
#                 Component Files - Upload and Data Checks                     #
################################################################################

# manifest file path
manifest_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/Manifests/ManifestsComplete")
# platemap file path
plate_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/PlateMaps/PlateMapsComplete")

# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

# nextclade file path
nc_fp <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

### output location for files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary")

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


mani_plate <- mani_plate %>% mutate(loc_code = case_when(received_source == "CDCIVY" ~ "IVY",
                                                         received_source == "CDCIVY5" ~ "IVY",
                                                         received_source == "CDCIVY6" ~ "IVY",
                                          received_source == "CDCIVY7" ~ "IVY",
                                          T ~ "UM"))


genbank <- read.csv(paste0(genbank_fp, "/sample_full_genbank_list.csv"), colClasses = "character")

mani_plate <- mani_plate %>% mutate(loc_code2 = case_when(received_source == "CDCIVY" ~  "IVY",
                                                          received_source == "CDCIVY5" ~ "IVY",
                                                          received_source == "CDCIVY6" ~ "IVY",
                                                  received_source == "CDCIVY7" ~ "IVY",
                                                  T ~ "UM"))

mani_plate_g2 <- merge(mani_plate, genbank, by.x = c("sample_id", "loc_code2"), by.y = c("sample_id", "loc_code2"), all.x = TRUE)
#mani_plate_pang_g2 <- merge(mani_plate_pang_g, gisaid, by.x = c("sample_id", "loc_code"), by.y = c("sample_id", "loc_code"), all.x = TRUE)


mppnc1 <- mani_plate_g2

nextclade <- read.csv(paste0(nextclade_fp, "/sample_full_nextclade_list.csv"), colClasses = "character")

mppnc2 <- merge(mppnc1, nextclade, by.x = c("sample_id", "PlateDate"), by.y = c("SampleID", "nextclade_runDate"), all.x = TRUE)

################################################################################
## add a column to number multiple sample_ids per subject_id

#if a failed run was processed through round2 of full_run_code.R
#then the line below will need to be uncommented to pull out the failed run data before running 
#gisaid_upload_file_creatin.R and the rest of this code will then need to be run
#fill in the PlateDat and PlateNumber of the run that you DON'T want the data from
#mppnc2 <- filter(mppnc2, PlateDate != "2022-06-06" | PlateNumber != "179")

mppnc2 <- mppnc2 %>% group_by(subject_id) %>% arrange(coll_date) %>% mutate(sample_per_subject = row_number())

# add a column for indicating if a particular subject_id has multiple samples
mppnc2 <- mppnc2 %>% group_by(subject_id) %>% mutate(multiSamples = case_when(max(sample_per_subject) > 1 ~ 1, 
                                                                              T ~ 0))
# filter out multiple samples
multiple_samples <- filter(mppnc2, multiSamples == 1)
# want to have an indicator marking samples that are 90 days from previous
multiple_samples <- multiple_samples %>% group_by(subject_id) %>% arrange(coll_date) %>% 
                      mutate(daysFromPrevious = as.numeric(as_date(coll_date) - lag(as_date(coll_date), default = NA)), 
                             ninetyDayFromPrevious = case_when(daysFromPrevious >= 90 ~ 1, 
                                                               T ~ 0), 
                             previousLineageDifferentThanCurrent = case_when(pangolin_lineage != lag(pangolin_lineage, default = NA) ~ 1, 
                                                                             T ~ 0), 
                             previousCladeDifferentThanCurrent = case_when(nextclade_clade != lag(nextclade_clade, default = NA) ~ 1, 
                                                                             T ~ 0))

mppnc2 <- merge(mppnc2, multiple_samples, all.x = TRUE)


mppnc2 <- mppnc2 %>% select(sample_id, subject_id, coll_date,                   
                                                       flag, received_source, SampleBarcode,               
                                                       PlateDate, PlatePlatform, PlateNumber,                 
                                                       genbank_SequenceID, genbank_Accession, genbank_SubmissionID,
                                                       received_date, position, SiteName,                    
                                                       subject_id_length, PlateName, PlatePosition,               
                                                       SampleSourceLocation, 
                                                       sample_per_subject, 
                                                       multiSamples, daysFromPrevious, ninetyDayFromPrevious)

mppnc2 <- mppnc2 %>% mutate(coll_date = case_when(grepl("/", coll_date) & substr(coll_date, nchar(coll_date) - 2, nchar(coll_date) - 2) != "/" ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                                  grepl("/", coll_date) & substr(coll_date, nchar(coll_date) - 2, nchar(coll_date) - 2) == "/" ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%y")), 
                                                      grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                                      T ~ NA_character_))


################################################################################

### add in data quality rule
#mppnc2 <- mppnc2 %>% mutate(data_quality_rule = case_when((pangolin_status %in% c("pass", "passed_qc")) & (nextclade_qcOverallStatus %in% c("good", "mediocre")) & (nextclade_completeness > 80) ~ "pass", 
#                                                          T ~ "not passed"))



################################################################################
## read in the clade.csv file and match to the data
clade_assign <- read.csv(paste0(starting_path, "SEQUENCING/hMPV/", plate_name, "_clades.csv"))


mppnc2 <- merge(mppnc2, clade_assign, by.x = "sample_id", by.y = "Sample", all.x = TRUE)

################################################################################
### negative control well warning

neg_control <- unique(filter(mppnc2, grepl("NC_", sample_id) & is.na(as.numeric(sample_id)))$sample_id)
neg_control2 <- unique(filter(mppnc2, grepl("NC-", sample_id) & is.na(as.numeric(sample_id)))$sample_id)
helas <- unique(filter(mppnc2, grepl("hela", tolower(sample_id)))$sample_id)

check_NCs <- filter(mppnc2, sample_id %in% neg_control | sample_id %in% helas | sample_id %in% neg_control2)

# We want to make sure with each plate that the three negative controls have ???10% of genome covered. 

check_NCs <- check_NCs %>% mutate(neg_control_warning = case_when(as.numeric(nextclade_completeness) >= 10 ~ 1,
                                                                  T ~ 0))

keep_NCs <- table(check_NCs$PlateName, check_NCs$neg_control_warning)

write.table(keep_NCs, paste0(outputLOC, "/ReportNotifications/negative_control_warnings.tsv"), sep = "\t")

################################################################################

write.csv(mppnc2, paste0(outputLOC, "/full_compiled_data.csv"), row.names = FALSE, na = "")
write.csv(mppnc2, paste0(outputLOC, "/secret/full_compiled_data.csv"), row.names = FALSE, na = "")

# a <- filter(mppnc2, SF456L_present == 1) %>% group_by(received_source, coll_date) %>% summarize(total = sum(as.numeric(SF456L_present)))
# 
# ggplot(a, aes(x = as_date(coll_date), y = total, fill = received_source)) + 
#   geom_col(stat = "identity") + 
#   scale_x_date(breaks = "1 month", date_labels = "%b '%y") +
#   theme_bw() +
#   labs(x = "Collection Date", 
#        y = "Count of S:F456L Detection", 
#        fill = "Source") + 
#   geom_text(aes(x = as_date("2022-09-24"), y = 1, label = "ASC - 007013444")) + 
#   geom_text(aes(x = as_date("2022-03-11"), y = 1, label = "CDCIVY4 - A93J40A2"))


