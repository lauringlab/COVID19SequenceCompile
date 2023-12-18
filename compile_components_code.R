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

# genbank file path
genbank_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/SequenceOutcomeComplete")


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

gisaid_secret <- filter(gisaid, grepl("RVTN", gisaid_strain) | grepl("VIEW", gisaid_strain))
#gisaid <- filter(gisaid, !grepl("RVTN", gisaid_strain))

# merge these onto mani_plate, always keeping everything from mani_plate
mani_plate_pang <- merge(mani_plate, pangolin, by.x = c("sample_id", "PlateDate"), by.y = c("SampleID", "pangolin_runDate"), all.x = TRUE)

if (nrow(mani_plate_pang) > nrow(mani_plate)){
  stop("Merging of Pangolin Data onto mainfest+plate maps = too many rows")
}

#mani_plate_pang <- filter(mani_plate_pang, received_source != "RVTN")
#mani_plate_pang_secret <- filter(mani_plate_pang, received_source == "RVTN")
### add column for time in days from plate to pangolin
#mani_plate_pang$PlateToPangolin_days <- difftime(mani_plate_pang$pangolin_runDate, mani_plate_pang$PlateDate, units = "days")

#table(mani_plate_pang$received_source)

mani_plate_pang <- mani_plate_pang %>% mutate(loc_code = case_when(received_source == "CDCIVY" ~ "IVY",
                                          received_source == "CDCIVY4" ~ "IVY",
                                          received_source == "CDCIVY5" ~ "IVY",
                                          received_source == "CDCIVY6" ~ "IVY",
                                          received_source == "RVTN" ~ "RVTN",
                                          received_source == "VIEW" ~ "VIEW",
                                          received_source == "IVYIC" ~ "IVYIC",
                                          T ~ "UM"))


mani_plate_pang_g <- merge(mani_plate_pang, gisaid, by.x = c("sample_id", "loc_code"), by.y = c("sample_id", "loc_code"), all.x = TRUE)
#mani_plate_pang_g_secret <- merge(mani_plate_pang_secret, gisaid_secret, by.x = c(" "), by.y = c("sample_id"), all.x = TRUE)

genbank <- read.csv(paste0(genbank_fp, "/sample_full_genbank_list.csv"), colClasses = "character")

genbank_secret <- filter(genbank, grepl("RVTN", genbank_SequenceID) | grepl("VIEW", genbank_SequenceID))

mani_plate_pang_g <- mani_plate_pang_g %>% mutate(loc_code2 = case_when(received_source == "CDCIVY" ~  "IVY",
                                                  received_source == "CDCIVY4" ~ "IVY",
                                                  received_source == "CDCIVY5" ~ "IVY",
                                                  received_source == "CDCIVY6" ~ "IVY",
                                                  received_source == "RVTN" ~ "RVTN",
                                                  received_source == "VIEW" ~ "VIEW",
                                                  received_source == "IVYIC" ~ "IVYIC",
                                                  received_source == "HFHS" ~ "MIS",
                                                  received_source == "ASC" ~ "MIS",
                                                  received_source == "ASJ" ~ "MIS",
                                                  received_source == "TRINITY" ~ "MIS",
                                                  received_source == "MDHHS" ~ "UM",
                                                  T ~ "UM"))

mani_plate_pang_g2 <- merge(mani_plate_pang_g, genbank, by.x = c("sample_id", "loc_code2"), by.y = c("sample_id", "loc_code2"), all.x = TRUE)
#mani_plate_pang_g2 <- merge(mani_plate_pang_g, gisaid, by.x = c("sample_id", "loc_code"), by.y = c("sample_id", "loc_code"), all.x = TRUE)


mppnc <- merge(mani_plate_pang_g2, nextclade, by.x = c("sample_id", "PlateDate"), by.y = c("SampleID", "nextclade_runDate"), all.x = TRUE)

if (nrow(mppnc) > nrow(mani_plate_pang_g )){
  stop("Merging of Nextclade Data onto mainfest+plate maps+pangolin+gisaid = too many rows")
}


### add column for time in days from plate to nextclade
#mppnc$PlateToNextclade_days <- difftime(mppnc$nextclade_runDate, mppnc$PlateDate, units = "days")

#### read in data from Emily's MHOME stuff
#mhome_in <- read.csv("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/10_transfer/MHome_HIVE/together.csv")
mhome_in <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/10_transfer/MHome_HIVE/together.csv"))
mhome_in$loc_code <- "UM"
mhome_in$loc_code2 <- "UM"
mhome_in$genbank_Accession <- ""
mhome_in$genbank_SequenceID <- ""
mhome_in$genbank_SubmissionID <- ""
mhome_in$SF456L_present <- ""
mhome_in <- mhome_in %>% select(colnames(mppnc))

mppnc <- rbind(mppnc, mhome_in)

################################################################################
# create indicator for if Plate to Nextclade or Plate to Pangolin is out of 
# expected range. this will help detect potential "wrong matches" for sample_ids 
# that are sent to us twice & re-run through process

# mppnc$IlluminaPangolin_OutOfRange <- ifelse(mppnc$PlatePlatform == "Illumina" & mppnc$PlateToPangolin_days > 8, 1, 0)
# mppnc$NanoporePangolin_OutOfRange <- ifelse(mppnc$PlatePlatform == "Nanopore" & mppnc$PlateToPangolin_days > 4, 1, 0)
# mppnc$IlluminaNextclade_OutOfRange <- ifelse(mppnc$PlatePlatform == "Illumina" & mppnc$PlateToNextclade_days > 8, 1, 0)
# mppnc$NanoporeNextclade_OutOfRange <- ifelse(mppnc$PlatePlatform == "Nanopore" & mppnc$PlateToNextclade_days > 4, 1, 0)
# 
# outofrangeset <- filter(mppnc, IlluminaPangolin_OutOfRange == 1 | NanoporePangolin_OutOfRange == 1 | IlluminaNextclade_OutOfRange == 1 | NanoporeNextclade_OutOfRange == 1)
# 
# outofrange_output <- filter(mppnc, sample_id %in% outofrangeset$sample_id)
# 
# write.csv(outofrange_output, paste0(outputLOC, "/ReportNotifications/out_of_range_alert.csv"), row.names = FALSE, na = "")

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

################################################################################
## apply logic for mis-matched pangolin/nextclade info
## necessary for instances where a sample was run on two different plates

# mppnc_look <- filter(mppnc2, sample_id %in% unique(filter(mppnc2, sample_per_subject > 1)$sample_id))
# mppnc_look$correct_matched <- ifelse(mppnc_look$PlateDate == mppnc_look$pangolin_runDate & mppnc_look$PlateDate == mppnc_look$nextclade_runDate, 1, 0)
# 
#a <- filter(mppnc2, sample_id == "10041282602")
# 
# mppnc2 <- mppnc2 %>% mutate(correct_matched = case_when(PlateDate == pangolin_runDate & PlateDate == nextclade_runDate ~ 1, 
#                                                         T ~ 0))
# 
# mppnc2 <- mppnc2 %>% group_by(sample_id) %>% mutate(count_platedates = length(unique(PlateDate)), 
#                                                     count_platedates2 = length(PlateDate),
#                                                     sum_matched = sum(correct_matched, na.rm = T))
# 
# mppnc2 <- mppnc2 %>% mutate(keeps = case_when(count_platedates == sum_matched ~ 1, 
#                                               T ~ 0))
# 
# mppnc2_outs <- filter(mppnc2, count_platedates2 %% 4 == 0 & sample_id != "")
# mppnc2_outs_keep <- filter(mppnc2_outs, correct_matched == 1)
# 
# goal <- nrow(mppnc2) - nrow(mppnc2_outs)
# mppnc2_t <- filter(mppnc2, count_platedates2 %% 4 != 0 | (count_platedates2 %% 4 == 0 & sample_id == ""))
# 
# if (nrow(mppnc2_t) != goal){
#   stop("Filter & Keep did not work properly")
# }
# 
# mppnc2 <- rbind(mppnc2_t, mppnc2_outs_keep) %>% select(sample_id, subject_id, coll_date,                   
#                                                        flag, received_source, SampleBarcode,               
#                                                        PlateDate, PlatePlatform, PlateNumber,                 
#                                                        pangolin_lineage, pangolin_probability, pangolin_status,             
#                                                        pangolin_note, nextclade_clade, nextclade_totalMissing,      
#                                                        nextclade_completeness, gisaid_strain, gisaid_epi_isl, 
#                                                        gisaid_clade, gisaid_pango_lineage, 
#                                                        received_date, position, SiteName,                    
#                                                        subject_id_length, PlateName, PlatePosition,               
#                                                        SampleSourceLocation, pangoLEARN_version, pangolin_conflict,           
#                                                        pango_version, pangolin_version, pangolin_runDate,            
#                                                        #PlateToPangolin_days, 
#                                                        nextclade_qcOverallScore, nextclade_qcOverallStatus,  
#                                                        nextclade_totalMutations, nextclade_totalNonACGTNs, nextclade_runDate,          
#                                                        #PlateToNextclade_days, IlluminaPangolin_OutOfRange, NanoporePangolin_OutOfRange, 
#                                                        #IlluminaNextclade_OutOfRange, NanoporeNextclade_OutOfRange, 
#                                                        sample_per_subject, 
#                                                        multiSamples, daysFromPrevious, ninetyDayFromPrevious, previousLineageDifferentThanCurrent, 
#                                                        previousCladeDifferentThanCurrent)


mppnc2 <- mppnc2 %>% select(sample_id, subject_id, coll_date,                   
                                                       flag, received_source, SampleBarcode,               
                                                       PlateDate, PlatePlatform, PlateNumber,                 
                                                       pangolin_lineage, pangolin_probability, pangolin_status,             
                                                       pangolin_note, nextclade_clade, nextclade_totalMissing,      
                                                       nextclade_completeness, gisaid_strain, gisaid_epi_isl, 
                                                       gisaid_clade, gisaid_pango_lineage, 
                                                       genbank_SequenceID, genbank_Accession, genbank_SubmissionID,
                                                       received_date, position, SiteName,                    
                                                       subject_id_length, PlateName, PlatePosition,               
                                                       SampleSourceLocation, pangoLEARN_version, pangolin_conflict,           
                                                       pango_version, pangolin_version, #pangolin_runDate,            
                                                       #PlateToPangolin_days, 
                                                       nextclade_qcOverallScore, nextclade_qcOverallStatus,  
                                                       nextclade_totalMutations, nextclade_totalNonACGTNs, SF456L_present, #nextclade_runDate,          
                                                       #PlateToNextclade_days, IlluminaPangolin_OutOfRange, NanoporePangolin_OutOfRange, 
                                                       #IlluminaNextclade_OutOfRange, NanoporeNextclade_OutOfRange, 
                                                       sample_per_subject, 
                                                       multiSamples, daysFromPrevious, ninetyDayFromPrevious, previousLineageDifferentThanCurrent, 
                                                       previousCladeDifferentThanCurrent)

mppnc2 <- mppnc2 %>% mutate(coll_date = case_when(grepl("/", coll_date) & substr(coll_date, nchar(coll_date) - 2, nchar(coll_date) - 2) != "/" ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                                  grepl("/", coll_date) & substr(coll_date, nchar(coll_date) - 2, nchar(coll_date) - 2) == "/" ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%y")), 
                                                      grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                                      T ~ NA_character_))

################################################################################

###
# pull in covid RVTN data
seq <- read.csv(paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"))

# only keep RVTN and VIEW
seq <- filter(seq, received_source %in% c("RVTN", "VIEW"))

# read in already assigned sequences
already_assigned <- read.csv(paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/RVTN/SampleID_Hide/assigned_rvtn_random.csv"))
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

write.csv(full_set_complete, paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/RVTN/SampleID_Hide/assigned_rvtn_random.csv"), row.names = FALSE, na = "")


# read in and attach RVTN and VIEW re-codes
rvtn_recodes <- read.csv(paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/RVTN/SampleID_Hide/assigned_rvtn_random.csv"), colClasses = "character")
rvtn_recodes <- rvtn_recodes %>% select(sample_id_lauring, sample_id)
rvtn_recodes <- filter(rvtn_recodes, sample_id != "")
#colnames(rvtn_recodes)
#colnames(mppnc2)

mppnc2 <- merge(mppnc2, rvtn_recodes, by = c("sample_id"), all.x = TRUE)

################################################################################

# add in RVTN gisaid
mppnc2_rvtn <- filter(mppnc2, grepl("RVTN", received_source))# received_source == "RVTN")
mppnc2_view <- filter(mppnc2, grepl("VIEW", received_source))
mppnc2 <- filter(mppnc2,!grepl("RVTN", received_source) & !grepl("VIEW", received_source)) #received_source != "RVTN")

mppnc2_rvtn <- mppnc2_rvtn %>% select(subject_id, sample_id, coll_date, flag,                               
                                      received_source, SampleBarcode,                     
                                      PlateDate, PlatePlatform,                      
                                      PlateNumber, pangolin_lineage,                  
                                      pangolin_probability, pangolin_status,                    
                                      pangolin_note, nextclade_clade,                    
                                      nextclade_totalMissing, nextclade_completeness,
                                      genbank_SequenceID, genbank_Accession, genbank_SubmissionID,
                                      received_date, position,                           
                                      SiteName, subject_id_length,                  
                                      PlateName, PlatePosition,                      
                                      SampleSourceLocation, pangoLEARN_version,                
                                      pangolin_conflict, pango_version,                     
                                      pangolin_version, 
                                      #pangolin_runDate,                   
                                      nextclade_qcOverallScore, nextclade_qcOverallStatus,          
                                      nextclade_totalMutations, nextclade_totalNonACGTNs, SF456L_present,          
                                      #nextclade_runDate, 
                                      sample_per_subject,                 
                                      multiSamples, daysFromPrevious,                  
                                      ninetyDayFromPrevious, previousLineageDifferentThanCurrent,
                                      previousCladeDifferentThanCurrent, sample_id_lauring)

mppnc2_view <- mppnc2_view %>% select(subject_id, sample_id, coll_date, flag,                               
                                      received_source, SampleBarcode,                     
                                      PlateDate, PlatePlatform,                      
                                      PlateNumber, pangolin_lineage,                  
                                      pangolin_probability, pangolin_status,                    
                                      pangolin_note, nextclade_clade,                    
                                      nextclade_totalMissing, nextclade_completeness,
                                      #genbank_SequenceID, genbank_Accession, genbank_SubmissionID,
                                      received_date, position,                           
                                      SiteName, subject_id_length,                  
                                      PlateName, PlatePosition,                      
                                      SampleSourceLocation, pangoLEARN_version,                
                                      pangolin_conflict, pango_version,                     
                                      pangolin_version, 
                                      #pangolin_runDate,                   
                                      nextclade_qcOverallScore, nextclade_qcOverallStatus,          
                                      nextclade_totalMutations, nextclade_totalNonACGTNs, SF456L_present,         
                                      #nextclade_runDate, 
                                      sample_per_subject,                 
                                      multiSamples, daysFromPrevious,                  
                                      ninetyDayFromPrevious, previousLineageDifferentThanCurrent,
                                      previousCladeDifferentThanCurrent, sample_id_lauring)

mppnc2_rvtn <- mppnc2_rvtn %>% mutate(loc_code = case_when(received_source == "CDCIVY" ~ "IVY",
                                                                   received_source == "CDCIVY4" ~ "IVY",
                                                                   received_source == "CDCIVY5" ~ "IVY",
                                                                   received_source == "CDCIVY6" ~ "IVY",
                                                                   received_source == "RVTN" ~ "RVTN",
                                                                   received_source == "VIEW" ~ "VIEW",
                                                                   received_source == "IVYIC" ~ "IVYIC",
                                                                   T ~ "UM"))

mppnc2_view <- mppnc2_view %>% mutate(loc_code = case_when(received_source == "CDCIVY" ~ "IVY",
                                                                   received_source == "CDCIVY4" ~ "IVY",
                                                                   received_source == "CDCIVY5" ~ "IVY",
                                                                   received_source == "CDCIVY6" ~ "IVY",
                                                                   received_source == "RVTN" ~ "RVTN",
                                                                   received_source == "VIEW" ~ "VIEW",
                                                                   received_source == "IVYIC" ~ "IVYIC",
                                                                   T ~ "UM"))

mppnc2_rvtn <- merge(mppnc2_rvtn, gisaid_secret, by.x = c("sample_id_lauring", "loc_code"), by.y = c("sample_id", "loc_code"), all.x = TRUE)
mppnc2_view <- merge(mppnc2_view, gisaid_secret, by.x = c("sample_id_lauring", "loc_code"), by.y = c("sample_id", "loc_code"), all.x = TRUE)

mppnc2_view <- merge(mppnc2_view, genbank_secret, by.x = c("sample_id_lauring", "loc_code"), by.y = c("sample_id", "loc_code2"), all.x = TRUE)

mppnc2_rvtn <- mppnc2_rvtn %>% select(subject_id, sample_id, coll_date, flag,                               
                                      received_source, SampleBarcode,                     
                                      PlateDate, PlatePlatform,                      
                                      PlateNumber, pangolin_lineage,                  
                                      pangolin_probability, pangolin_status,                    
                                      pangolin_note, nextclade_clade,                    
                                      nextclade_totalMissing, nextclade_completeness,             
                                      gisaid_strain, gisaid_epi_isl,                     
                                      gisaid_clade, gisaid_pango_lineage,  
                                      genbank_SequenceID, genbank_Accession, genbank_SubmissionID,
                                      received_date, position,                           
                                      SiteName, subject_id_length,                  
                                      PlateName, PlatePosition,                      
                                      SampleSourceLocation, pangoLEARN_version,                
                                      pangolin_conflict, pango_version,                     
                                      pangolin_version, 
                                      #pangolin_runDate,                   
                                      nextclade_qcOverallScore, nextclade_qcOverallStatus,          
                                      nextclade_totalMutations, nextclade_totalNonACGTNs, SF456L_present,         
                                      #nextclade_runDate, 
                                      sample_per_subject,                 
                                      multiSamples, daysFromPrevious,                  
                                      ninetyDayFromPrevious, previousLineageDifferentThanCurrent,
                                      previousCladeDifferentThanCurrent, sample_id_lauring)

mppnc2_view <- mppnc2_view %>% select(subject_id, sample_id, coll_date, flag,                               
                                      received_source, SampleBarcode,                     
                                      PlateDate, PlatePlatform,                      
                                      PlateNumber, pangolin_lineage,                  
                                      pangolin_probability, pangolin_status,                    
                                      pangolin_note, nextclade_clade,                    
                                      nextclade_totalMissing, nextclade_completeness,             
                                      gisaid_strain, gisaid_epi_isl,                     
                                      gisaid_clade, gisaid_pango_lineage,
                                      genbank_SequenceID, genbank_Accession, genbank_SubmissionID,
                                      received_date, position,                           
                                      SiteName, subject_id_length,                  
                                      PlateName, PlatePosition,                      
                                      SampleSourceLocation, pangoLEARN_version,                
                                      pangolin_conflict, pango_version,                     
                                      pangolin_version, 
                                      #pangolin_runDate,                   
                                      nextclade_qcOverallScore, nextclade_qcOverallStatus,          
                                      nextclade_totalMutations, nextclade_totalNonACGTNs, SF456L_present,         
                                      #nextclade_runDate, 
                                      sample_per_subject,                 
                                      multiSamples, daysFromPrevious,                  
                                      ninetyDayFromPrevious, previousLineageDifferentThanCurrent,
                                      previousCladeDifferentThanCurrent, sample_id_lauring)


mppnc2 <- rbind(mppnc2, mppnc2_rvtn)
mppnc2 <- rbind(mppnc2, mppnc2_view)

rm(mppnc2_rvtn)
rm(mppnc2_view)

################################################################################

### add in data quality rule
mppnc2 <- mppnc2 %>% mutate(data_quality_rule = case_when((pangolin_status %in% c("pass", "passed_qc")) & (nextclade_qcOverallStatus %in% c("good", "mediocre")) & (nextclade_completeness > 80) ~ "pass", 
                                                          T ~ "not passed"))

################################################################################

# get full pangolin file
#/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate
full_pangolin_new <- list.files(path = paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate"), pattern = "lineage_report*")
fpn <- read.csv(paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/pangolin/CompleteFastaUpToDate/", full_pangolin_new))

# only select the sample id and lineage call
fpn <- fpn %>% select(taxon, lineage)
colnames(fpn) <- c("sample_id", "newest_pangolin_lineage")

# pull the date portion out and attach that
date_bit <- substr(full_pangolin_new, 16, 23)
fpn$newest_pangolin_date <- date_bit

fpn <- fpn %>% group_by(sample_id) %>% mutate(count = length(sample_id)) %>% distinct()
#fpn <- filter(fpn, count == 1)
fpn <- fpn %>% select(sample_id, newest_pangolin_lineage, newest_pangolin_date)

### remove out any negative controls, etc.
fpn <- filter(fpn, !grepl("NC_", sample_id) & !grepl("HeLa", sample_id) & !grepl("NC-", sample_id))

# merge that data onto full set
mppnc2 <- merge(mppnc2, fpn, by = c("sample_id"), all.x = TRUE, all.y = FALSE)

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


