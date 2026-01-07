### Code to generate RedCap upload file for RSVA and RSVB - RIGHT
### Created: 06/03/2024
### Author: Leigh Papalambros
#################################################

library(tidyverse)
library(lubridate)
library(dplyr)


# read in full compiled rsv_a file

rsv_a <- read.csv("/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character", stringsAsFactors = FALSE)

right1_a <- filter(rsv_a, received_source == "RIGHT")

right1_a$subtype <- "A"
#right1_a_date_try <- data.frame()
#right1_a_date_try <- as.Date(right1_a$coll_date, tryFormats = c("%y-%m-%d", "%Y-%m-%d"))

#right1_a$coll_date <- as.POSIXct(right1_a$coll_date, format = c("%y-%m-%d", "%Y-%m-%d"))

rsv_b <- read.csv("/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character", stringsAsFactors = FALSE)

right1_b <- filter(rsv_b, received_source == "RIGHT")

right1_b$subtype <- "B"

right1_ab <- rbind(right1_a, right1_b)
str(right1_ab$coll_date)

#right1_ab$coll_date <- as.Date(as.POSIXct(right1_ab$coll_date, tryFormats = c("%y-%m-%d", "%Y-%m-%d", "%m-%d-%y")))

samp_dup <- right1_ab %>% group_by(subject_id) %>% summarize(count = n ()) %>% filter(count > 1)
print(samp_dup)
dups_list <- merge(samp_dup, right1_ab, by = "subject_id", all.x = TRUE)

#write.csv(dups_list, "~/Desktop/right_rsv_dups.csv", row.names = FALSE)

#non_dups_for_upload <- anti_join(right1_ab, samp_dup, by = "subject_id")

# replacing the blank nextclade cells with NA per RIGHT request

right1_ab <- right1_ab %>% mutate(sample_id = case_when(sample_id == "308309N01" ~ "308402N01",
                                                         sample_id == "308310N03" ~ "308403N03",
                                                         sample_id == "308311N03" ~ "308404N03",
                                                        sample_id == "308310N02" ~ "308403N02",
                                                         T ~ sample_id))

## editing the dupe samples to add "_2" and "_3" to the samples 
right1_ab <- right1_ab %>% mutate(subject_id = case_when(subject_id == "307903" & sample_id == "307903N07" ~ "307903_2",
                                                    subject_id == "307902" & sample_id == "307902N06" ~ "307902_2",
                                                    subject_id == "308001" & sample_id == "308001N02" ~ "308001_2",
                                                    subject_id == "308001" & sample_id == "308001N03" ~ "308001_3",
                                                    subject_id == "307901" & sample_id == "307901N05" ~ "307901_2",
                                                    subject_id == "307803" & sample_id == "307803N04" ~ "307803_2",
                                                    subject_id == "307505" & sample_id == "307505N02" ~ "307505_2",
                                                    subject_id == "307503" & sample_id == "307503N02" ~ "307503_2",
                                                    subject_id == "307401" & sample_id == "307401N02" ~ "307401_2",
                                                    subject_id == "221501" & sample_id == "221501N10" ~ "221501_2",
                                                    subject_id == "220802" & sample_id == "220802N04" ~ "220802_2",
                                                    subject_id == "220801" & sample_id == "220801N04" ~ "220801_2",
                                                    subject_id == "220602" & sample_id == "220602N07" ~ "220602_2",
                                                    subject_id == "219601" & sample_id == "219601N02" ~ "219601_2",
                                                    subject_id == "219301" & sample_id == "219301N08" ~ "219301_2",
                                                    subject_id == "218905" & sample_id == "218905N05" ~ "218905_2",
                                                    subject_id == "218901" & sample_id == "218901N10" ~ "218901_2",
                                                    subject_id == "218203" & sample_id == "218203N09" ~ "218203_2",
                                                    subject_id == "214001" & sample_id == "214001N03" ~ "214001_2",
                                                    subject_id == "213801" & sample_id == "213801N05" ~ "213801_2",
                                                    subject_id == "213304" & sample_id == "213304N02" ~ "213304_2",
                                                    subject_id == "213301" & sample_id == "213301N10" ~ "213301_2",
                                                    subject_id == "212702" & sample_id == "212702N03" ~ "212702_2",
                                                    subject_id == "211801" & sample_id == "211801N09" ~ "211801_2",
                                                    subject_id == "211001" & sample_id == "211001N03" ~ "211001_2",
                                                    subject_id == "210201" & sample_id == "210201N05" ~ "210201_2",
                                                    subject_id == "126001" & sample_id == "126001N03" ~ "126001_2",
                                                    subject_id == "125102" & sample_id == "125102N05" ~ "125102_2",
                                                    subject_id == "123301" & sample_id == "123301N02" ~ "123301_2",
                                                    subject_id == "106501" & sample_id == "106501N02" ~ "106501_2",
                                                    subject_id == "308301" & sample_id == "308301N03" ~ "308301_2",
                                                    subject_id == "308310" & sample_id == "308403N03" ~ "308310_2",
                                                    subject_id == "308802" & sample_id == "308802N02" ~ "308802_2",
                                                    subject_id == "309001" & sample_id == "309001N04" ~ "309001_2",
                                                    subject_id == "223801" & sample_id == "223801N06" ~ "223801_2",
                                                    subject_id == "223802" & sample_id == "223802N07" ~ "223802_2",
                                                    subject_id == "123401" & sample_id == "123401N05" ~ "123401_2",
                                                    subject_id == "123501" & sample_id == "123501N04" ~ "123501_2",
                                                    subject_id == "123503" & sample_id == "123503N09" ~ "123503_2",
                                                    subject_id == "124605" & sample_id == "124605N03" ~ "124605_2",
                                                    subject_id == "222801" & sample_id == "222801N04" ~ "222801_2",
                                                    T ~ subject_id))


right1_ab <- right1_ab %>% mutate_at(c('nextclade_clade', 'nextclade_Gclade', 
                                  'nextclade_totalMutations', 'nextclade_totalNonACGTNs', 
                                  'nextclade_totalMissing', 'nextclade_completeness', 'nextclade_qcOverallStatus'), ~na_if(., ''))
#print(df)

#sapply(right1_ab_out, class)

# three samples need sample_id changes, they were hand written and the cdc wants their ID's



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



# double checking if there are any more dupes not corrected
samp_dup_2 <- right1_ab_out %>% group_by(subject_id) %>% summarize(count = n ()) %>% filter(count > 1)
print(samp_dup_2)
#dups_list_2 <- merge(samp_dup_2, right1_ab_out, by = "subject_id", all.x = TRUE)


f_out <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/RIGHT_uploads/"


write.csv(right1_ab_out, paste0(f_out, "right_rsv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")




