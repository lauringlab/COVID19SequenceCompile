#### Influenza CDC IVY Upload

library(tidyverse)
library(lubridate)

#### Read in influenza full file
flu_file <- read.csv(paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/", 
                            "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                            "full_compiled_data.csv"), colClasses = c("character"))

flu_fileB <- read.csv(paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/", 
                            "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/", 
                            "full_compiled_data.csv"), colClasses = c("character"))

rv_flu <- filter(flu_file, grepl("RVTN", received_source) & !grepl("Missing Date", flag))
rv_fluB <- filter(flu_fileB, grepl("RVTN", received_source) & !grepl("Missing Date", flag))

colnames(rv_flu)

#### segment ids in as blanks currently
# cdc_flu <- cdc_flu %>% mutate(PB2.Segment_ID = "", 
#                               PB1.Segment_ID = "", 
#                               PA.Segment_ID = "", 
#                               NP.Segment_ID = "", 
#                               NA.Segment_ID = "", 
#                               MP.Segment_ID = "", 
#                               NS.Segment_ID = "")


rv_flu <- rv_flu %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                              received_date, SampleBarcode, PlateName, nextclade_HA_clade, 
                              nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                              nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                              nextclade_HA_type, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                              NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, Isolate_Name)

rv_fluB <- rv_fluB %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                            received_date, SampleBarcode, PlateName, nextclade_HA_clade, 
                            nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                            nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                            nextclade_HA_type, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                            NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, Isolate_Name)

rv_flu <- rbind(rv_flu, rv_fluB)

colnames(rv_flu) <- c("sample_id", "subject_id", "coll_date_flu", "flag_flu", 
                       "received_source_flu", "received_date_flu", "sample_barcode_flu", 
                       "plate_name_flu", "nextclade_ha_clade_flu", "nextclade_ha_qcoverallscore_flu", 
                       "nextclade_ha_qcoverallstatus_flu", "nextclade_ha_totalmutations_flu", 
                       "nextclade_ha_totalnonacgtns_flu", "nextclade_ha_type_flu", "isolate_id_flu", "pb2_segment_id_flu", 
                       "pb1_segment_id_flu", "pa_segment_id_flu", "ha_segment_id_flu", "np_segment_id_flu", 
                       "na_segment_id_flu", "mp_segment_id_flu", "ns_segment_id_flu", "isolate_name_flu")


# table(rv_flu$plate_name_flu)
# rv_flu_out <- filter(rv_flu, plate_name_flu == "20230315_IAV_Illumina_Run_49")
# rv_flu_out[, c(9:24)] <- ""
# rv_flu_out$flag_flu <- "Removed - Failed Negative Control Well Check"
# rv_flu <- filter(rv_flu, plate_name_flu != "20230315_IAV_Illumina_Run_49")
# flu_all <- rbind(rv_flu, rv_flu_out)


today_date <- gsub("-", "", Sys.Date())
write.csv(rv_flu, paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/", 
                          "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                          "RVTN_uploads/rvtn_flu_", today_date, ".csv"), row.names = FALSE, na = "")
