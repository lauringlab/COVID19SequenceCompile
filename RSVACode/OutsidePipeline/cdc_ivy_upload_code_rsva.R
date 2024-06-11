### Code to generate RedCap upload file for RSVA and RSVB - IVY


library(tidyverse)
library(lubridate)


# read in full compiled rsv_a file

rsv_a <- read.csv("/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

rsv_a4 <- filter(rsv_a, received_source == "CDCIVY4")
rsv_a5 <- filter(rsv_a, received_source == "CDCIVY5")
rsv_a6 <- filter(rsv_a, received_source == "CDCIVY6")

rsv_a4$subtype <- "A"
rsv_a5$subtype <- "A"
rsv_a6$subtype <- "A"

rsv_b <- read.csv("/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

rsv_b4 <- filter(rsv_b, received_source == "CDCIVY4")
rsv_b5 <- filter(rsv_b, received_source == "CDCIVY5")
rsv_b6 <- filter(rsv_b, received_source == "CDCIVY6")

rsv_b4$subtype <- "B"
rsv_b5$subtype <- "B"
rsv_b6$subtype <- "B"

rsv_ab4 <- rbind(rsv_a4, rsv_b4)
rsv_ab5 <- rbind(rsv_a5, rsv_b5)
rsv_ab6 <- rbind(rsv_a6, rsv_b6)

rsv_ab_out4 <- rsv_ab4 %>% select(subject_id, sample_id, coll_date, flag, 
                              received_source, received_date, position, 
                              PlateName, PlateDate, PlatePosition, SampleSourceLocation, 
                              subtype, genbank_SubmissionID, genbank_Accession, genbank_SequenceID,
                              nextclade_clade, nextclade_Gclade, 
                              nextclade_totalMutations, nextclade_totalNonACGTNs, 
                              nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)

rsv_ab_out5 <- rsv_ab5 %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, position, 
                                  PlateName, PlatePosition, SampleSourceLocation, 
                                  subtype, epi_isl,
                                  nextclade_clade, nextclade_Gclade, 
                                  nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                  nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)
sapply(rsv_ab_out5, class)

rsv_ab_out6 <- rsv_ab6 %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, position, 
                                  PlateName, PlateDate, PlatePosition, SampleSourceLocation, 
                                  subtype, genbank_SubmissionID, genbank_Accession, genbank_SequenceID,
                                  nextclade_clade, nextclade_Gclade, 
                                  nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                  nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)


colnames(rsv_ab_out4) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                         "received_source_rsv", "received_date_rsv", "position_rsv", 
                         "plate_name_rsv", "plate_date_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                         "subtype_rsv", "genbank_submissionid_rsv", "genbank_accessionid_rsv", "genbank_sequenceid_rsv",
                         "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                         "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                         "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")

colnames(rsv_ab_out5) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                           "received_source_rsv", "received_date_rsv", "position_rsv", 
                           "plate_name_rsv","plate_position_rsv", "sample_source_location_rsv", 
                           "subtype_rsv", "gisaid_epi_isl_rsv",
                           "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                           "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                           "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")

colnames(rsv_ab_out6) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                           "received_source_rsv", "received_date_rsv", "position_rsv", 
                           "plate_name_rsv", "plate_date_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                           "subtype_rsv", "genbank_submissionid_rsv", "genbank_accessionid_rsv", "genbank_sequenceid_rsv",
                           "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                           "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                           "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")


f_out <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/IVY_uploads/"


#write.csv(rsv_ab_out4, paste0(f_out, "ivy4_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")
#write.csv(rsv_ab_out5, paste0(f_out, "ivy5_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")
write.csv(rsv_ab_out6, paste0(f_out, "ivy6_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")




