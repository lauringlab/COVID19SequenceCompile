### Code to generate RedCap upload file for RSVA - IVY

library(tidyverse)
library(lubridate)


# read in full compiled rsv_a file

rsv_a <- read.csv("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

rsv_a <- filter(rsv_a, received_source == "CDCIVY5")

rsv_a$subtype <- "A"

rsv_a_out <- rsv_a %>% select(subject_id, sample_id, coll_date, flag, 
                              received_source, received_date, position, 
                              PlateName, PlatePosition, SampleSourceLocation, 
                              epi_isl, nextclade_clade, nextclade_Gclade, 
                              subtype, nextclade_totalMutations, nextclade_totalNonACGTNs, 
                              nextclade_totalMissing, nextclade_completeness)


colnames(rsv_a_out) <- c("subject_id", "sample_id", "coll_date_rsv", "flag_rsv", 
                         "received_source_rsv", "received_date_rsv", "position_rsv", 
                         "plate_name_rsv", "plate_position_rsv", "sample_source_location_rsv", 
                         "gisaid_epi_isl_rsv", "nextclade_clade_rsv", "nextclade_gclade_rsv", 
                         "subtype_rsv", "nextclade_totalmutations_rsv", "nextclade_totalnonacgtns_rsv", 
                         "nextclade_totalmissing_rsv", "nextclade_completeness_rsv")


f_out <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/IVY_uploads/"

write.csv(rsv_a_out, paste0(f_out, "ivy5_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")

