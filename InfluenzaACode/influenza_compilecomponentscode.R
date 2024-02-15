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
gisaid_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/gisaid")
# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")

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


plate_map <- filter(plate_map, !grepl("IBV", SampleSourceLocation))

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

#if (nrow(mani_plate) != plate_map_ids){
#  stop("There are more or less rows in our manifest + plate combination than there were date/sample id combinations in the original plate map file")
#}

#missings <- filter(mani_plate, is.na(subject_id))
#write.csv(missings, "C:/Users/juliegil/Documents/UofM_Work/Lauring_Lab/check_miss_subjects.csv", na = "", row.names = FALSE)

# then, read in pangolin, gisaid, and next clade files
nextclade <- read.csv(paste0(nc_fp, "/sample_full_nextclade_list.csv"), colClasses = "character")


### gisaid
### read every .csv in that folder
g_files <- list.files(gisaid_fp, pattern = "*.csv")

gisaid <- data.frame()

for (i in g_files){
  gisaid_in <- read.csv(paste0(gisaid_fp, "/", i), colClasses = "character")
  gisaid <- rbind(gisaid, gisaid_in)
}

colnames(gisaid)[1] <- "Isolate_Id"

gisaid_secret <- filter(gisaid, grepl("RVTN", Isolate_Name))
gisaid <- filter(gisaid, !grepl("RVTN", Isolate_Name))

gisaid <- gisaid %>% select(Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                            NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id, 
                            P3.Segment_Id, Isolate_Name)
gisaid_secret <- gisaid_secret %>% select(Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                            NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id, 
                            P3.Segment_Id, Isolate_Name)

gisaid <- separate(data = gisaid, col = Isolate_Name, remove = FALSE, sep = "/", into = c("type", "place", "id", "year"))
gisaid_secret <- separate(data = gisaid_secret, col = Isolate_Name, remove = FALSE, sep = "/", into = c("type", "place", "id", "year"))

gisaid$id <- gsub("UOM", "", gisaid$id)
gisaid$id <- gsub("IVY", "", gisaid$id)
gisaid$id <- gsub("HFHS", "", gisaid$id)
gisaid$id <- gsub("RVTN", "", gisaid$id)
gisaid$id <- gsub("MISAPPHIRE", "", gisaid$id)
gisaid_secret$id <- gsub("RVTN", "", gisaid_secret$id)

gisaid <- gisaid %>% select(id, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                            NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id, 
                            P3.Segment_Id, Isolate_Name)
colnames(gisaid) <- c("sample_id", "Isolate_Id", "PB2.Segment_Id", "PB1.Segment_Id", "PA.Segment_Id", "HA.Segment_Id", 
                      "NP.Segment_Id", "NA.Segment_Id", "MP.Segment_Id", "NS.Segment_Id", "HE.Segment_Id", 
                      "P3.Segment_Id", "Isolate_Name")

gisaid_secret <- gisaid_secret %>% select(id, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                            NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id, 
                            P3.Segment_Id, Isolate_Name)
colnames(gisaid_secret) <- c("sample_id", "Isolate_Id", "PB2.Segment_Id", "PB1.Segment_Id", "PA.Segment_Id", "HA.Segment_Id", 
                      "NP.Segment_Id", "NA.Segment_Id", "MP.Segment_Id", "NS.Segment_Id", "HE.Segment_Id", 
                      "P3.Segment_Id", "Isolate_Name")


mani_plate_g <- merge(mani_plate, gisaid, by.x = c("sample_id"), by.y = c("sample_id"), all.x = TRUE)

################################################################################
### genbank

genbank <- read.csv(paste0(genbank_fp, "/sample_full_genbank_list.csv"), colClasses = "character")

colnames(genbank)[1] <- "genbank_SubmissionID"

genbank_secret <- filter(genbank, grepl("RVTN", loc_code2))
genbank <- filter(genbank, !grepl("RVTN", loc_code2))

mani_plate_g2 <- merge(mani_plate_g, genbank, by.x = c("sample_id"), by.y = c("sample_id"), all.x = TRUE)

################################################################################

mppnc <- merge(mani_plate_g2, nextclade, by.x = c("sample_id"), by.y = c("SampleID"), all.x = TRUE)

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
### rvtn recode set-up

###
# pull in flu RVTN data
seq <- read.csv(paste0(starting_path, "/SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"))

# only keep RVTN
seq <- filter(seq, received_source == "RVTN")

# read in already assigned sequences
already_assigned <- read.csv(paste0(starting_path, "/SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/RVTN/SampleID_Hide/assigned_rvtn_random.csv"))
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

full_set2 <- cbind(not_assigned, seq3) ## this contains all newly assigned rvtn stuff, plus all the unassigned ids

full_set_complete <- rbind(filter(already_assigned, !is.na(subject_id)), full_set2)

write.csv(full_set_complete, paste0(starting_path, "/SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/RVTN/SampleID_Hide/assigned_rvtn_random.csv"), row.names = FALSE, na = "")


# read in and attach RVTN re-codes
rvtn_recodes <- read.csv(paste0(starting_path, "/SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/RVTN/SampleID_Hide/assigned_rvtn_random.csv"), colClasses = "character")
rvtn_recodes <- rvtn_recodes %>% select(sample_id_lauring, sample_id, subject_id)
#colnames(rvtn_recodes)
#colnames(mppnc2)

mppnc2 <- merge(mppnc, rvtn_recodes, by = c("subject_id", "sample_id"), all.x = TRUE)

################################################################################

# add in RVTN gisaid
mppnc2_rvtn <- filter(mppnc2, grepl("RVTN", received_source))# received_source == "RVTN")
mppnc2 <- filter(mppnc2,!grepl("RVTN", received_source)) #received_source != "RVTN")
 
mppnc2_rvtn <- mppnc2_rvtn %>% select(sample_id, subject_id, coll_date,                   
                                      flag, received_source, received_date, SampleBarcode,
                                      PlateDate, PlatePlatform, PlateNumber,
                                      genbank_SubmissionID, loc_code2, genbank_HA, genbank_HAH1,               
                                      genbank_HAH3, genbank_MP, genbank_NA, genbank_NAN1,              
                                      genbank_NAN2, genbank_NP, genbank_NS,                 
                                      genbank_PA, genbank_PB1, genbank_PB2,
                                      nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,
                                      nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus,
                                      nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                                      nextclade_HA_runDate, nextclade_HA_type,
                                      subject_id_length, position, PlateName, PlatePosition,
                                      SampleSourceLocation, PlateToNextclade_days,
                                      sample_per_subject, sample_id_lauring)
 
mppnc2_rvtn <- merge(mppnc2_rvtn, gisaid_secret, by.x = c("sample_id_lauring"), by.y = c("sample_id"), all.x = TRUE)

#colnames(mppnc2_rvtn)
colnames(mppnc2)
#colnames(gisaid_secret)

mppnc2_rvtn <- mppnc2_rvtn %>% select(sample_id, subject_id, coll_date,
                                      flag, received_source, received_date, SampleBarcode,
                                      PlateDate, PlatePlatform, PlateNumber,
                                      genbank_SubmissionID, loc_code2, genbank_HA, genbank_HAH1,               
                                      genbank_HAH3, genbank_MP, genbank_NA, genbank_NAN1,              
                                      genbank_NAN2, genbank_NP, genbank_NS,                 
                                      genbank_PA, genbank_PB1, genbank_PB2,
                                      nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,
                                      nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus,
                                      nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                                      nextclade_HA_runDate, nextclade_HA_type,
                                      Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id,
                                      NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id,
                                      P3.Segment_Id, Isolate_Name,
                                      subject_id_length, position, PlateName, PlatePosition,
                                      SampleSourceLocation, PlateToNextclade_days,
                                      sample_per_subject, sample_id_lauring)

mppnc2 <- rbind(mppnc2, mppnc2_rvtn)

rm(mppnc2_rvtn)


################################################################################

# add in RVTN genbank
mppnc2_rvtn <- filter(mppnc2, grepl("RVTN", received_source))# received_source == "RVTN")
mppnc2 <- filter(mppnc2,!grepl("RVTN", received_source)) #received_source != "RVTN")

mppnc2_rvtn <- mppnc2_rvtn %>% select(sample_id, subject_id, coll_date,                   
                                      flag, received_source, received_date, SampleBarcode,
                                      PlateDate, PlatePlatform, PlateNumber,
                                      nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,
                                      nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus,
                                      nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                                      nextclade_HA_runDate, nextclade_HA_type,
                                      Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id,
                                      NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id,
                                      P3.Segment_Id, Isolate_Name,
                                      subject_id_length, position, PlateName, PlatePosition,
                                      SampleSourceLocation, PlateToNextclade_days,
                                      sample_per_subject, sample_id_lauring)

mppnc2_rvtn <- merge(mppnc2_rvtn, genbank_secret, by.x = c("sample_id_lauring"), by.y = c("sample_id"), all.x = TRUE)

#colnames(mppnc2_rvtn)

mppnc2_rvtn <- mppnc2_rvtn %>% select(sample_id, subject_id, coll_date,
                                      flag, received_source, received_date, SampleBarcode,
                                      PlateDate, PlatePlatform, PlateNumber,
                                      genbank_SubmissionID, loc_code2, genbank_HA, genbank_HAH1,               
                                      genbank_HAH3, genbank_MP, genbank_NA, genbank_NAN1,              
                                      genbank_NAN2, genbank_NP, genbank_NS,                 
                                      genbank_PA, genbank_PB1, genbank_PB2,
                                      nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,
                                      nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus,
                                      nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                                      nextclade_HA_runDate, nextclade_HA_type,
                                      Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id,
                                      NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id,
                                      P3.Segment_Id, Isolate_Name,
                                      subject_id_length, position, PlateName, PlatePosition,
                                      SampleSourceLocation, PlateToNextclade_days,
                                      sample_per_subject, sample_id_lauring)

mppnc2 <- rbind(mppnc2, mppnc2_rvtn)

rm(mppnc2_rvtn)

################################################################################
## creating "StrainName" for genbank submissions
## this is equivalent to "Isolate_Name" from old gisiad submission ways

#getting state information for IVY
#if (any(grepl("IVY", mppnc2$received_source))){
  
#  mppnc2$source_state_code <- substr(mppnc2$subject_id, 3, 4)
#  cdcivy_manifest_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/CDCIVY")
#  cdc_sites <- read.csv(paste0(cdcivy_manifest_fp, "/Keys/CDC_SiteCodebook.csv"), colClasses = "character") %>% separate(SiteCode, into = c("site", "state"), sep = "_")
#  site_bit <- cdc_sites %>% select(Number, state)
  # add leading zeros to site
#  site_bit <- site_bit %>% mutate(Number = case_when(nchar(Number) == 1 ~ paste0("0", Number), 
#                                                     T ~ Number))
  #colnames(site_bit) <- c("Number", "state")
#  mppnc2 <- merge(mppnc2, site_bit, by.x = c("source_state_code"), by.y = c("Number"), all.x = TRUE)
#  mppnc2$state <- ifelse(!grepl("IVY", mppnc2$received_source), "", mppnc2$state)
  
#} else {
#  state <- ""
#}

## need to add in other options of Michigan for state if not IVY above

## creating StrainName variable for the various sources
## A/State/genbank_id/Year
#mppnc2 <- mppnc2 %>% mutate(StrainName = case_when(received_source == "CDCIVY6" ~ paste0("A/", state, "/IVY-", genbank_id, "/", substr(coll_date, 1, 4)),
#                                                                                       received_source == "CDCIVY4" ~ paste0("SARS-CoV-2/human/USA/", "IVY-", sample_id, "/", substr(coll_date, 1, 4)),


################################################################################
## apply logic for mis-matched pangolin/nextclade info
## necessary for instances where a sample was run on two different plates
## unsure how this logic will work for flu, since not everything gets a nextclade entry

#colnames(mppnc)


mppnc3 <- mppnc2 %>% select(sample_id, subject_id, coll_date,                   
                           flag, received_source, received_date, SampleBarcode,               
                           PlateDate, PlatePlatform, PlateNumber,
                           genbank_SubmissionID, loc_code2, genbank_HAH1,               
                           genbank_HAH3, genbank_MP, genbank_NAN1,              
                           genbank_NAN2, genbank_NP, genbank_NS,                 
                           genbank_PA, genbank_PB1, genbank_PB2,
                           nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,      
                           nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                           nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                           nextclade_HA_runDate, nextclade_HA_type, 
                           Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                           NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, HE.Segment_Id, 
                           P3.Segment_Id, Isolate_Name,              
                           subject_id_length, position, PlateName, PlatePosition,               
                           SampleSourceLocation, PlateToNextclade_days, 
                           sample_per_subject, sample_id_lauring)

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

write.csv(mppnc3, paste0(outputLOC, "/full_compiled_data.csv"), row.names = FALSE, na = "")
write.csv(mppnc3, paste0(outputLOC, "/secret/full_compiled_data.csv"), row.names = FALSE, na = "")
