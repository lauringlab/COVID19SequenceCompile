################################################################################
#       Creation of Cumulative Dataset for RSVB Genetic Sampling          #
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
manifest_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/ManifestsComplete")
# platemap file path
plate_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/PlateMaps/PlateMapsComplete")
# nextclade file path
nc_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
# gisaid file path
gisaid_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/SequenceOutcomes/gisaid")
# previous 2021 file path
#prev_2021 <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PreviousLists")
# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")
### output location for files, all together
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary")

################################################################################


# first, read in manifest file
manifest <- read.csv(paste0(manifest_fp, "/sample_full_manifest_list.csv"), colClasses = "character")

# read in plate map file
plate_map <- read.csv(paste0(plate_fp, "/sample_full_plate_list.csv"), colClasses = "character")
plate_map <- filter(plate_map, SampleID != "" & !is.na(SampleID)) # remove any rows where sample ID is missing on the plate map
#plate_map$SampleID <- as.character(trimws(plate_map$SampleID))


#plate_map <- filter(plate_map, !grepl("IBV", SampleSourceLocation))

# store unique number of sample ids 
plate_map_ids <- nrow(plate_map %>% group_by(SampleID, SampleSourceDate) %>% summarize(count = length(SampleID)))


################################################################################
# Warning for if plate map date and manifest date are DIFFERENT

manifest_options <- filter(manifest, received_date != "" & !is.na(received_date)) %>% select(sample_id, subject_id, received_date, received_source)
platemap_options <- filter(plate_map, SampleSourceDate != "" & !is.na(SampleSourceDate)) %>% select(SampleID, SampleSourceDate, PlateNumber, SampleSourceLocation)

compare_options <- merge(manifest_options, platemap_options, by.x = c("sample_id"), by.y = c("SampleID"))

compare_options$different <- ifelse(compare_options$received_date != compare_options$SampleSourceDate, 1, 0)

compare_options <- filter(compare_options, as_date(received_date) >= as_date("2021-07-01") & different == 1)

if (nrow(compare_options) > 0){
  print(compare_options)
  stop("Mismatched received dates between manifest and platemap.")
}

################################################################################

manifest$sample_id <- trimws(manifest$sample_id)
plate_map$SampleID <- trimws(plate_map$SampleID)
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

################################################################################

# then, read in gisaid and next clade files
nextclade <- read.csv(paste0(nc_fp, "/sample_full_nextclade_list.csv"), colClasses = "character")


### gisaid
### read every .csv in that folder
# g_files <- list.files(gisaid_fp, pattern = "*.csv")
# 
# gisaid <- data.frame()
# 
# for (i in g_files){
#   gisaid_in <- read.csv(paste0(gisaid_fp, "/", i), colClasses = "character")
#   gisaid <- rbind(gisaid, gisaid_in)
# }
# 
# colnames(gisaid)[1] <- "Isolate_Id"

gisaid <- read.csv(paste0(gisaid_fp, "/gisaid_back_info.csv"), colClasses = "character")
gisaid <- gisaid %>% select(sample_id, epi_isl, strain_name)



mani_plate_g <- merge(mani_plate, gisaid, by.x = c("sample_id"), by.y = c("sample_id"), all.x = TRUE)
mppnc <- merge(mani_plate_g, nextclade, by.x = c("sample_id"), by.y = c("SampleID"), all.x = TRUE)

genbank <- read.csv(paste0(genbank_fp, "/sample_full_genbank_list.csv"), colClasses = "character")

mppnc <- merge(mppnc, genbank, by.x = c("sample_id"), by.y = c("sample_id"), all.x = TRUE)
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
## adding in the lauring_lab_id recode info

genbank_secret <- filter(genbank, grepl("RIGHT", genbank_SequenceID))

# pull in covid RVTN, VIEW, and RIGHT data
seq <- read.csv(paste0(starting_path, "/SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"))

# only keep RIGHT
seq <- filter(seq, received_source %in% c("RIGHT"))

# read in already assigned sequences
already_assigned <- read.csv(paste0(starting_path, "/SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/RIGHT/SampleID_Hide/assigned_right_random.csv"))
already_assigned <- already_assigned %>% select(sample_id_lauring, sample_id, subject_id)

### only keep items in seq that are NOT already assigned
seq2 <- filter(seq, !sample_id %in% unique(already_assigned$sample_id))

# filter out already_assigned so only non-assigned lauring labels are present
not_assigned <- filter(already_assigned, is.na(subject_id)) %>% select(sample_id_lauring)

# pull out sample & subject id, add to full_set
seq3 <- seq2 %>% select(subject_id, sample_id)
fillup <- data.frame(rep(NA, nrow(not_assigned)-nrow(seq3)), rep(NA, nrow(not_assigned)-nrow(seq3)))
colnames(fillup) <- colnames(seq3)
seq3 <- rbind(seq3, fillup)

full_set2 <- cbind(not_assigned, seq3) ## this contains all newly assigned right stuff, plus all the unassigned ids

full_set_complete <- rbind(filter(already_assigned, !is.na(subject_id)), full_set2)

write.csv(full_set_complete, paste0(starting_path, "/SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/RIGHT/SampleID_Hide/assigned_rvtn_random.csv"), row.names = FALSE, na = "")


# read in and attach RIGHT re-codes
right_recodes <- read.csv(paste0(starting_path, "/SEQUENCING/RSV_B/4_SequenceSampleMetadata/Manifests/RIGHT/SampleID_Hide/assigned_rvtn_random.csv"), colClasses = "character")
right_recodes <- right_recodes %>% select(sample_id_lauring, sample_id)
right_recodes <- filter(right_recodes, sample_id != "")
#colnames(right_recodes)


mppnc <- merge(mppnc, right_recodes, by = c("sample_id"), all.x = TRUE)




################################################################################

mppnc3 <- mppnc 
# %>% select(sample_id, subject_id, coll_date,                   
#                            flag, received_source, received_date, SampleBarcode,               
#                            PlateDate, PlatePlatform, PlateNumber,
#                            nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,      
#                            nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
#                            nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
#                            nextclade_HA_runDate, nextclade_HA_type, 
#                            Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
#                            NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id, 
#                            P3.Segment_Id, Isolate_Name,              
#                            subject_id_length, position, PlateName, PlatePosition,               
#                            SampleSourceLocation, PlateToNextclade_days, 
#                            sample_per_subject, sample_id_lauring)

## removing duplicate IVY6 samples that were sequenced with the wrong primers


################################################################################
### negative control well warning

neg_control <- unique(filter(mppnc3, grepl("NC_", sample_id) & is.na(as.numeric(sample_id)))$sample_id)
neg_control2 <- unique(filter(mppnc3, grepl("NC-", sample_id) & is.na(as.numeric(sample_id)))$sample_id)
helas <- unique(filter(mppnc3, grepl("hela", tolower(sample_id)))$sample_id)

check_NCs <- filter(mppnc3, sample_id %in% neg_control | sample_id %in% helas | sample_id %in% neg_control2)

# We want to make sure with each plate that the three negative controls have ???10% of genome covered. 

check_NCs <- check_NCs %>% mutate(neg_control_warning = case_when(as.numeric(nextclade_completeness) >= 10 ~ 1,
                                                                  T ~ 0))

keep_NCs <- table(check_NCs$PlateName, check_NCs$neg_control_warning)

write.table(keep_NCs, paste0(outputLOC, "/ReportNotifications/negative_control_warnings.tsv"), sep = "\t")
################################################################################

write.csv(mppnc3, paste0(outputLOC, "/full_compiled_data.csv"), row.names = FALSE, na = "")
write.csv(mppnc3, paste0(outputLOC, "/secret/full_compiled_data.csv"), row.names = FALSE, na = "")
