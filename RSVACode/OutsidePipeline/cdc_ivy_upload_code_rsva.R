### Code to generate RedCap upload file for RSVA - IVY

library(tidyverse)
library(lubridate)


# read in full compiled rsv_a file

rsv_a <- read.csv("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

rsv_a4 <- filter(rsv_a, received_source == "CDCIVY4")
rsv_a5 <- filter(rsv_a, received_source == "CDCIVY5")

rsv_a4$subtype <- "A"
rsv_a5$subtype <- "A"

rsv_b <- read.csv("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

rsv_b4 <- filter(rsv_b, received_source == "CDCIVY4")
rsv_b5 <- filter(rsv_b, received_source == "CDCIVY5")

rsv_b4$subtype <- "B"
rsv_b5$subtype <- "B"

rsv_a4 <- rbind(rsv_a4, rsv_b4)
rsv_a5 <- rbind(rsv_a5, rsv_b5)

rsv_a_out4 <- rsv_a4 %>% select(subject_id, sample_id, coll_date, flag, 
                              received_source, received_date, position, 
                              PlateName, PlatePosition, SampleSourceLocation, 
                              epi_isl, nextclade_clade, nextclade_Gclade, 
                              subtype, nextclade_totalMutations, nextclade_totalNonACGTNs, 
                              nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)

rsv_a_out5 <- rsv_a5 %>% select(subject_id, sample_id, coll_date, flag, 
                                received_source, received_date, position, 
                                PlateName, PlatePosition, SampleSourceLocation, 
                                epi_isl, nextclade_clade, nextclade_Gclade, 
                                subtype, nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallStatus)


colnames(rsv_a_out4) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                         "received_source_rsv", "received_date_rsv", "position_rsv", 
                         "plate_name_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                         "gisaid_epi_isl_rsv", "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                         "subtype_rsv", "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                         "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")

colnames(rsv_a_out5) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                          "received_source_rsv", "received_date_rsv", "position_rsv", 
                          "plate_name_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                          "gisaid_epi_isl_rsv", "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                          "subtype_rsv", "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                          "nextclade_totalmissing_rsv", "nextclade_completeness_rsv", "nextclade_overall_qc_rsv")


f_out <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/IVY_uploads/"


write.csv(rsv_a_out4, paste0(f_out, "ivy4_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")
write.csv(rsv_a_out5, paste0(f_out, "ivy5_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")

