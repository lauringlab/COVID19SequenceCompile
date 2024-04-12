#### Influenza CDC IVY Upload
library(tidyverse)
library(lubridate)

#### Read in influenza full file
flu_file <- read.csv(paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/", 
                            "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                            "full_compiled_data.csv"), colClasses = c("character"))

flu_fileB <- read.csv(paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/", 
                            "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/", 
                            "full_compiled_data.csv"), colClasses = c("character"))
#flu_file <- filter(flu_file, PlateName != "20230315_IAV_Illumina_Run_49")


cdc_flu <- filter(flu_file, grepl("IVY", received_source) & !grepl("Missing Date", flag))
cdc_fluB <- filter(flu_fileB, grepl("IVY", received_source) & !grepl("Missing Date", flag))
table(cdc_flu$received_source)
table(cdc_fluB$received_source)
colnames(cdc_flu)

#### segment ids in as blanks currently
# cdc_flu <- cdc_flu %>% mutate(PB2.Segment_ID = "", 
#                               PB1.Segment_ID = "", 
#                               PA.Segment_ID = "", 
#                               NP.Segment_ID = "", 
#                               NA.Segment_ID = "", 
#                               MP.Segment_ID = "", 
#                               NS.Segment_ID = "")

#cdc_fluB <- cdc_fluB %>% mutate(Isolate_Id = "",
#                               PB2.Segment_ID = "", 
#                               PB1.Segment_ID = "", 
#                               PA.Segment_ID = "", 
#                               HA.Segment_Id = "",
#                               NP.Segment_ID = "", 
#                               NA.Segment_ID = "", 
#                              MP.Segment_ID = "", 
#                              NS.Segment_ID = "",
#                              Isolate_Name = "")

cdc_flu_old <- cdc_flu %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                              received_date, SampleBarcode, PlateName, nextclade_HA_clade, 
                              nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                              nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                              nextclade_HA_type, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                              NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, Isolate_Name)

cdc_fluB_old <- cdc_fluB %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                              received_date, SampleBarcode, PlateName, nextclade_HA_clade, 
                              nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                              nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                              nextclade_HA_type, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                              NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, Isolate_Name)

cdc_flu_old <- rbind(cdc_flu_old, cdc_fluB_old)


cdc_flu_new <- cdc_flu %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                                  received_date, SampleBarcode, PlateName, PlateDate, PlatePlatform, PlateNumber,
                                  PlatePosition, SampleSourceLocation, nextclade_HA_clade, nextclade_HA_type,
                                  nextclade_HA_completeness, nextclade_HA_totalMissing, nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                                  nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                                  genbank_SubmissionID, genbank_HAH1, genbank_HAH3,
                                  genbank_MP, genbank_NAN1, genbank_NAN2, genbank_NP, genbank_NS,
                                  genbank_PA, genbank_PB1, genbank_PB2)

cdc_fluB$genbank_HAH1 <- ""
cdc_fluB$genbank_HAH3 <- ""
cdc_fluB$genbank_NAN1 <- ""
cdc_fluB$genbank_NAN2 <- ""

cdc_fluB_new <- cdc_fluB %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                                    received_date, SampleBarcode, PlateName, PlateDate, PlatePlatform, PlateNumber,
                                    PlatePosition, SampleSourceLocation, nextclade_HA_clade, nextclade_HA_type,
                                    nextclade_HA_completeness, nextclade_HA_totalMissing, nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                                    nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                                    genbank_SubmissionID, genbank_HAH1, genbank_HAH3,
                                    genbank_MP, genbank_NAN1, genbank_NAN2, genbank_NP, genbank_NS,
                                    genbank_PA, genbank_PB1, genbank_PB2)

cdc_flu_new_full <- rbind(cdc_flu_new, cdc_fluB_new)

# cdc_flu <- cdc_flu %>% mutate(site_code = substr(subject_id, 3, 4))
# 
# a <- filter(cdc_flu, site_code %in% c("23", "17"))
# a %>% group_by(received_source, site_code, substr(PlateName, 10, 12)) %>% summarize(count = length(unique(sample_id)))

colnames(cdc_flu_old) <- c("sample_id", "study_id", "coll_date_flu", "flag_flu", 
                       "received_source_flu", "received_date_flu", "sample_barcode_flu", 
                       "plate_name_flu", "nextclade_ha_clade_flu", "nextclade_ha_qcoverallscore_flu", 
                       "nextclade_ha_qcoverallstatus_flu", "nextclade_ha_totalmutations_flu", 
                       "nextclade_ha_totalnonacgtns_flu", "nextclade_ha_type_flu", "isolate_id_flu", "pb2_segment_id_flu", 
                       "pb1_segment_id_flu", "pa_segment_id_flu", "ha_segment_id_flu", "np_segment_id_flu", 
                       "na_segment_id_flu", "mp_segment_id_flu", "ns_segment_id_flu", "isolate_name_flu")

colnames(cdc_flu_new_full) <- c("sample_id", "subject_id", "coll_date_flu", "flag_flu", 
                           "received_source_flu", "received_date_flu", "sample_barcode_flu", 
                           "plate_name_flu", "plate_date_flu", "plate_platform_flu", "plate_number_flu",
                           "plate_position_flu", "sample_source_location_flu", "nextclade_ha_clade_flu", 
                           "nextclade_ha_type_flu", "nextclade_ha_completeness_flu", "nextclade_ha_totalmissing_flu",
                           "nextclade_ha_qcoverallscore_flu", "nextclade_ha_qcoverallstatus_flu", "nextclade_ha_totalmutations_flu", 
                           "nextclade_ha_totalnonacgtns_flu", "genbank_submissionid_flu", "genbank_hah1_flu",
                           "genbank_hah3_flu", "genbank_mp_flu", "genbank_nan1_flu", "genbank_nan2_flu", "genbank_np_flu",
                           "genbank_ns_flu", "genbank_pa_flu", "genbank_pb1_flu", "genbank_pb2_flu")

cdc_flu_old <- cdc_flu_old %>% mutate(study_id = case_when(sample_id == "G43Q59U6" ~ "23240105-01", 
                                                   T ~ study_id))

today_date <- gsub("-", "", Sys.Date())

#write.csv(filter(cdc_flu_old, received_source_flu == "CDCIVY4"), paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/", 
#                          "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
#                          "IVY_uploads/cdc_ivy_flu_", today_date, ".csv"), row.names = FALSE, na = "")


flu5 <- filter(cdc_flu_old, received_source_flu == "CDCIVY5")
# flu5_out <- filter(flu5, plate_name_flu %in% c("20230315_IAV_Illumina_Run_49", "20230316_IAV_Illumina_Run_51"))
# flu5_out[, c(9:24)] <- ""
# flu5_out$flag_flu <- "Removed - Failed Negative Control Well Check"
# flu5 <- filter(flu5, !plate_name_flu %in% c("20230315_IAV_Illumina_Run_49", "20230316_IAV_Illumina_Run_51"))
# flu_all <- rbind(flu5, flu5_out)

#write.csv(flu5, paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/", 
#                                                                    "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
#                                                                    "IVY_uploads/cdc_ivy5_flu_", today_date, ".csv"), row.names = FALSE, na = "")

flu6 <- filter(cdc_flu_new_full, received_source_flu == "CDCIVY6")

write.csv(flu6, paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/", 
                       "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                       "IVY_uploads/cdc_ivy6_flu_", today_date, ".csv"), row.names = FALSE, na = "")

