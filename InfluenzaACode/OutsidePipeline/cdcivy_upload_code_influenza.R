#### Influenza CDC IVY Upload
library(tidyverse)
library(lubridate)

#### Read in influenza full file
flu_file <- read.csv(paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/", 
                            "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                            "full_compiled_data.csv"), colClasses = c("character"))

cdc_flu <- filter(flu_file, grepl("IVY", received_source) & !grepl("Missing Date", flag))
table(cdc_flu$received_source)
colnames(cdc_flu)

#### segment ids in as blanks currently
# cdc_flu <- cdc_flu %>% mutate(PB2.Segment_ID = "", 
#                               PB1.Segment_ID = "", 
#                               PA.Segment_ID = "", 
#                               NP.Segment_ID = "", 
#                               NA.Segment_ID = "", 
#                               MP.Segment_ID = "", 
#                               NS.Segment_ID = "")


cdc_flu <- cdc_flu %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                              received_date, SampleBarcode, PlateName, nextclade_HA_clade, 
                              nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, 
                              nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs, 
                              nextclade_HA_type, Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, 
                              NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, NS.Segment_Id, Isolate_Name)

colnames(cdc_flu) <- c("sample_id", "study_id", "coll_date_flu", "flag_flu", 
                       "received_source_flu", "received_date_flu", "sample_barcode_flu", 
                       "plate_name_flu", "nextclade_ha_clade_flu", "nextclade_ha_qcoverallscore_flu", 
                       "nextclade_ha_qcoverallstatus_flu", "nextclade_ha_totalmutations_flu", 
                       "nextclade_ha_totalnonacgtns_flu", "nextclade_ha_type_flu", "isolate_id_flu", "pb2_segment_id_flu", 
                       "pb1_segment_id_flu", "pa_segment_id_flu", "ha_segment_id_flu", "np_segment_id_flu", 
                       "na_segment_id_flu", "mp_segment_id_flu", "ns_segment_id_flu", "isolate_name_flu")

today_date <- gsub("-", "", Sys.Date())

write.csv(filter(cdc_flu, received_source_flu == "CDCIVY4"), paste0("/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/", 
                          "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                          "IVY_uploads/cdc_ivy_flu_", today_date, ".csv"), row.names = FALSE, na = "")

write.csv(filter(cdc_flu, received_source_flu == "CDCIVY5"), paste0("/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/", 
                                                                    "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/", 
                                                                    "IVY_uploads/cdc_ivy5_flu_", today_date, ".csv"), row.names = FALSE, na = "")
