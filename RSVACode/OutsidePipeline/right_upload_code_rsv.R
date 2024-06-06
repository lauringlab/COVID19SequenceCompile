### Code to generate RedCap upload file for RSVA and RSVB - RIGHT
### Created: 06/03/2024
### Author: Leigh Papalambros
#################################################

library(tidyverse)
library(lubridate)


# read in full compiled rsv_a file

rsv_a <- read.csv("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

right1_a <- filter(rsv_a, received_source == "RIGHT")

right1_a$subtype <- "A"

rsv_b <- read.csv("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

right1_b <- filter(rsv_b, received_source == "RIGHT")

right1_b$subtype <- "B"

right1_ab <- rbind(right1_a, right1_b)

#sapply(right1_ab_out, class)

right1_ab_out <- right1_ab %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, position, 
                                  PlateName, PlateDate, PlatePosition, SampleSourceLocation, 
                                  subtype, genbank_SubmissionID, genbank_Accession, genbank_SequenceID,
                                  nextclade_clade, nextclade_Gclade, 
                                  nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                  nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)



colnames(right1_ab_out) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                           "received_source_rsv", "received_date_rsv", "position_rsv", 
                           "plate_name_rsv", "plate_date_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                           "subtype_rsv", "genbank_submissionid_rsv", "genbank_accessionid_rsv", "genbank_sequenceid_rsv",
                           "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                           "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                           "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")


f_out <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/RIGHT_uploads/"


write.csv(right1_ab_out, paste0(f_out, "right_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")




