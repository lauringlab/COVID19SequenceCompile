### Code to generate RedCap upload file for hMPV - IVY


library(tidyverse)
library(lubridate)


# read in full compiled rsv_a file

hMPV <- read.csv("/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

hmpv5 <- filter(hMPV, received_source == "CDCIVY5")
hmpv6 <- filter(hMPV, received_source == "CDCIVY6")
hpmv7 <- filter(hMPV, received_source == "CDCIVY7")


hmpv_out5 <- hmpv5 %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, position, 
                                  PlateName, PlatePosition, SampleSourceLocation, 
                                  subtype, epi_isl,
                                  nextclade_clade, nextclade_Gclade, 
                                  nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                  nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)
#sapply(rsv_ab_out5, class)

hmpv_out6 <- hmpvb6 %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, position, 
                                  PlateName, PlateDate, PlatePosition, SampleSourceLocation, 
                                  subtype, genbank_SubmissionID, genbank_Accession, genbank_SequenceID,
                                  nextclade_clade, nextclade_Gclade, 
                                  nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                  nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)

hmpv_out7 <- hmpv7 %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, position, 
                                  PlateName, PlateDate, PlatePosition, SampleSourceLocation, 
                                  subtype, genbank_SubmissionID, genbank_Accession, genbank_SequenceID,
                                  nextclade_clade, nextclade_Gclade, 
                                  nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                  nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)


colnames(hmpv_out5) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                           "received_source_rsv", "received_date_rsv", "position_rsv", 
                           "plate_name_rsv","plate_position_rsv", "sample_source_location_rsv", 
                           "subtype_rsv", "gisaid_epi_isl_rsv",
                           "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                           "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                           "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")

colnames(hmpv_out6) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                           "received_source_rsv", "received_date_rsv", "position_rsv", 
                           "plate_name_rsv", "plate_date_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                           "subtype_rsv", "genbank_submissionid_rsv", "genbank_accessionid_rsv", "genbank_sequenceid_rsv",
                           "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                           "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                           "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")

colnames(hmpv_out7) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                           "received_source_rsv", "received_date_rsv", "position_rsv", 
                           "plate_name_rsv", "plate_date_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                           "subtype_rsv", "genbank_submissionid_rsv", "genbank_accessionid_rsv", "genbank_sequenceid_rsv",
                           "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                           "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                           "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")


f_out <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary/IVY_uploads/"


write.csv(rsv_ab_out5, paste0(f_out, "ivy5_hmpv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")
write.csv(rsv_ab_out6, paste0(f_out, "ivy6_hmpv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")
write.csv(rsv_ab_out7, paste0(f_out, "ivy7_hmpv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")




