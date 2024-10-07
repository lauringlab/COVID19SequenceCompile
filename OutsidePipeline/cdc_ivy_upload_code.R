################################################################################
#       Creation of CDC IVY Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 07/20/2023                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################
checking_wd <- getwd()

if (grepl("leighbak", checking_wd)){
  starting_path <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/"
  outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/")

}else if (grepl("chbl", checking_wd)){
  starting_path <- "/Users/chbl/University of Michigan Dropbox/MED-LauringLab/"
  outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/")
} else {
  
  print("User not recognized.")
  
}


seq_list <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")
#seq_list_o <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################

# seq_list <- seq_list %>% select(sample_id, subject_id, coll_date,                    
#                                 flag, received_source, SiteName, SampleBarcode,                
#                                 PlateDate, PlatePlatform, PlateNumber,                 
#                                 pangolin_lineage, pangolin_probability, pangolin_status,             
#                                 pangolin_note, nextclade_clade, nextclade_totalMissing,      
#                                 nextclade_completeness, gisaid_strain, gisaid_epi_isl,              
#                                 received_date, position,          
#                                 PlateName, PlatePosition, SampleSourceLocation,        
#                                 pangoLEARN_version, pangolin_conflict, pango_version,               
#                                 pangolin_version, pangolin_runDate, #PlateToPangolin_days,        
#                                 nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,    
#                                 nextclade_totalNonACGTNs, nextclade_runDate)#, PlateToNextclade_days)

### for ivy 5
# seq_list <- seq_list %>% select(sample_id, subject_id, coll_date,                    
#                                 flag, received_source, SiteName, SampleBarcode,                
#                                 PlateDate, PlatePlatform, PlateNumber,                 
#                                 pangolin_lineage, pangolin_probability, pangolin_status,             
#                                 pangolin_note, nextclade_clade, nextclade_totalMissing,      
#                                 nextclade_completeness, gisaid_strain, gisaid_epi_isl,              
#                                 received_date, position,          
#                                 PlateName, PlatePosition, SampleSourceLocation,        
#                                 pangoLEARN_version, pangolin_conflict, pango_version,               
#                                 pangolin_version, pangolin_runDate,         
#                                 nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,    
#                                 nextclade_totalNonACGTNs, nextclade_runDate, data_quality_rule, newest_pangolin_lineage, newest_pangolin_date)

# seq_list_o <- seq_list_o %>% select(sample_id, subject_id, coll_date,
#                                flag, received_source, SiteName, SampleBarcode,
#                                PlateDate, PlatePlatform, PlateNumber,
#                                pangolin_lineage, pangolin_probability, pangolin_status,
#                                pangolin_note, nextclade_clade, nextclade_totalMissing,
#                                nextclade_completeness, gisaid_strain, gisaid_epi_isl,
# received_date, position,
# PlateName, PlatePosition, SampleSourceLocation,
# pangoLEARN_version, pangolin_conflict, pango_version,
# pangolin_version, pangolin_runDate, PlateToPangolin_days,
# nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,
# nextclade_totalNonACGTNs, nextclade_runDate, PlateToNextclade_days)


seq_list <- filter(seq_list, (received_source == "CDCIVY" | received_source == "CDCIVY4" | received_source == "CDCIVY5" | received_source == "CDCIVY6" | received_source == "CDCIVY7") & !grepl("Missing Date", flag))


## check for CDC IVY 4 samples (start with 22, ivy 3 == 21)
# if (any(substr(seq_list$subject_id, 1, 2) == 22)){
#   print("IVY 4 Samples Present, Need to separate for upload")
# }

if (length(unique(seq_list$sample_id)) != nrow(seq_list)){
  stop("Duplicate sample IDs - handle accordingly")
}

## remove study withdraws
seq_list <- filter(seq_list, flag != "Withdrawn from study")
seq_list <- filter(seq_list, !grepl("Withdrawn from study", flag))


################################################################################
### fix date formatting
seq_list <- seq_list %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

# seq_list_o <- seq_list_o %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")),
#                                                       grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")),
#                                                       T ~ NA_character_))


################################################################################
# 6/18/2021 - no longer need this - will be doing full overwrite upload going forward
# After the first upload, we'll need to keep track of what has already been uploaded
# so we'll read in the full list, then read in the previous upload list, and only keep 
# rows that are not in the previous upload list(s) to write out and upload the next time

# ivy_file_list <- list.files(pattern = "*.csv", path = paste0(outputLOC, "ARCHIVE/"))
# 
# ivy_redcap <- data.frame()
# for (i in ivy_file_list){
#   one <- read.csv(paste0(outputLOC, "ARCHIVE/", i), colClasses = "character")
#   ivy_redcap <- rbind(ivy_redcap, one)
# }
# 
# # select only rows from seqlist that are not in ivy_redcap
# combo <- anti_join(seq_list, ivy_redcap)

#length(unique(combo$sample_id))
#combo$flag <- ifelse(grepl("REDO", combo$SampleSourceLocation), "Re-run sample from batch #1", "")
################################################################################

colnames(seq_list) <- tolower(colnames(seq_list))

# add leading zero to month
if (nchar(month(Sys.Date())) == 1){
  m <- paste0("0", month(Sys.Date()))
} else {
  m <- month(Sys.Date())
}
# add leading zero to day
if (nchar(day(Sys.Date())) == 1){
  d <- paste0("0", day(Sys.Date()))
} else {
  d <- day(Sys.Date())
}

today <- paste0(year(Sys.Date()), m, d)


#ivy3 <- filter(seq_list, substr(subject_id, 1, 2) == 21)
ivy4 <- filter(seq_list, substr(subject_id, 1, 2) == 22)
ivy5 <- filter(seq_list, substr(subject_id, 1, 2) == 23)
ivy6 <- filter(seq_list, substr(subject_id, 1, 2) == 24)
#ivy7 <- filter(seq_list, substr(subject_id, 1, 2) == 24)


# change subject_id to study_id
ivy4 <- ivy4 %>% select(sample_id, subject_id, coll_date,                    
                                flag, received_source, sitename, samplebarcode,                
                                platedate, plateplatform, platenumber, 
                                pangolin_lineage, pangolin_status, pangolin_note,
                        nextclade_clade, nextclade_totalmissing, nextclade_completeness, 
                        gisaid_strain, gisaid_epi_isl, received_date, position, platename,
                        plateposition, samplesourcelocation, pangolearn_version,
                        pango_version, pangolin_version, nextclade_qcoverallscore, nextclade_qcoverallstatus, 
                                nextclade_totalmutations, nextclade_totalnonacgtns)

#SC2_NP_330
# ivy4_out <- filter(ivy4, platename %in% c("20220817_SC2_Illumina_Run_63"))
# ivy4_out[, c(11:30)] <- ""
# ivy4_out$flag <- "Removed - Failed Negative Control Well Check"
# ivy4 <- filter(ivy4, !platename %in% c("20220817_SC2_Illumina_Run_63"))
# ivy4 <- rbind(ivy4, ivy4_out)




names(ivy4)[names(ivy4) == 'subject_id'] <- 'study_id'


ivy5$nextclade_sf456l_present <- ivy5$sf456l_present

# change subject_id to study_id
ivy5 <- ivy5 %>% select(sample_id, subject_id, coll_date,                    
                        flag, received_source, sitename, samplebarcode,                
                        platedate, plateplatform, platenumber, 
                        pangolin_lineage, pangolin_status, pangolin_note,
                        nextclade_clade, nextclade_totalmissing, nextclade_completeness, 
                        gisaid_strain, gisaid_epi_isl, received_date, position, platename,
                        plateposition, samplesourcelocation, pangolearn_version,
                        pango_version, pangolin_version, nextclade_qcoverallscore, nextclade_qcoverallstatus, 
                        nextclade_totalmutations, nextclade_totalnonacgtns, 
                        data_quality_rule, newest_pangolin_lineage, newest_pangolin_date, nextclade_sf456l_present)

#SC2_NP_333, SC2_NP_334, SC2_NP_335

# ivy5 <- ivy5 %>% mutate(newest_pangolin_lineage = case_when(platename %in% c("20230720_SC2_Nanopore_Run_333", "20230720_SC2_Nanopore_Run_334", "20230720_SC2_Nanopore_Run_335") ~ "", 
#                                                             T ~ newest_pangolin_lineage), 
#                         newest_pangolin_date = case_when(platename %in% c("20230720_SC2_Nanopore_Run_333", "20230720_SC2_Nanopore_Run_334", "20230720_SC2_Nanopore_Run_335") ~ "", 
#                                                             T ~ newest_pangolin_date))

#table(ivy5$newest_pangolin_lineage == "XBB.1.5")
# ivy5_out <- filter(ivy5, platename %in% c("20230124_SC2_Illumina_Run_78", "20230330_SC2_Illumina_Run_95", "20230523_SC2_Illumina_Run_106"))
# ivy5_out[, c(11:33)] <- ""
# ivy5_out$flag <- "Removed - Failed Negative Control Well Check"
# ivy5 <- filter(ivy5, !platename %in% c("20230124_SC2_Illumina_Run_78", "20230330_SC2_Illumina_Run_95", "20230523_SC2_Illumina_Run_106"))
# ivy5 <- rbind(ivy5, ivy5_out)

#names(ivy5)[names(ivy5) == 'subject_id'] <- 'study_id'

ivy6 <- ivy6 %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                        received_date, sitename, samplebarcode, platename, platedate, 
                        plateplatform, platenumber, plateposition, samplesourcelocation, 
                        pangolin_lineage, pangolin_status, pangolin_note, 
                        pangolin_version, pangolin_conflict, nextclade_clade,
                        nextclade_totalmissing, nextclade_completeness, nextclade_qcoverallscore,
                        nextclade_qcoverallstatus, nextclade_totalmutations, 
                        nextclade_totalnonacgtns, genbank_sequenceid, genbank_accession, 
                        genbank_submissionid, data_quality_rule, newest_pangolin_lineage, 
                        newest_pangolin_date, sf456l_present)

colnames(ivy6) <- c("sample_id", "subject_id", "coll_date", "flag", "received_source", 
                    "received_date", "sitename", "samplebarcode", "platename", "platedate", 
                    "plateplatform", "platenumber", "plateposition", "samplesourcelocation", 
                    "pangolin_lineage", "pangolin_status", "pangolin_note", 
                    "pangolin_version", "pangolin_conflict", "nextclade_clade",
                    "nextclade_totalmissing", "nextclade_completeness", "nextclade_qcoverallscore",
                    "nextclade_qcoverallstatus", "nextclade_totalsubstitutions", 
                    "nextclade_totalnonacgtns", "genbank_sequenceid", "genbank_accession", 
                    "genbank_submissionid", "data_quality_rule", "newest_pangolin_lineage", 
                    "newest_pangolin_date", "nextclade_sf456l_present")


ivy7 <- ivy7 %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                        received_date, sitename, samplebarcode, platename, platedate, 
                        plateplatform, platenumber, plateposition, samplesourcelocation, 
                        pangolin_lineage, pangolin_status, pangolin_note, 
                        pangolin_version, pangolin_conflict, nextclade_clade,
                        nextclade_totalmissing, nextclade_completeness, nextclade_qcoverallscore,
                        nextclade_qcoverallstatus, nextclade_totalmutations, 
                        nextclade_totalnonacgtns, genbank_sequenceid, genbank_accession, 
                        genbank_submissionid, data_quality_rule, newest_pangolin_lineage, 
                        newest_pangolin_date, sf456l_present)

colnames(ivy7) <- c("sample_id", "subject_id", "coll_date", "flag", "received_source", 
                    "received_date", "sitename", "samplebarcode", "platename", "platedate", 
                    "plateplatform", "platenumber", "plateposition", "samplesourcelocation", 
                    "pangolin_lineage", "pangolin_status", "pangolin_note", 
                    "pangolin_version", "pangolin_conflict", "nextclade_clade",
                    "nextclade_totalmissing", "nextclade_completeness", "nextclade_qcoverallscore",
                    "nextclade_qcoverallstatus", "nextclade_totalsubstitutions", 
                    "nextclade_totalnonacgtns", "genbank_sequenceid", "genbank_accession", 
                    "genbank_submissionid", "data_quality_rule", "newest_pangolin_lineage", 
                    "newest_pangolin_date", "nextclade_sf456l_present")

#seq_list <- filter(seq_list, platenumber <= 49)
#write.csv(ivy3, paste0(outputLOC, "cdc_ivy3_", today, ".csv"), row.names = FALSE, na = "")
#write.csv(ivy4, paste0(outputLOC, "cdc_ivy4_", today, ".csv"), row.names = FALSE, na = "")
#write.csv(ivy5, paste0(outputLOC, "cdc_ivy5_", today, ".csv"), row.names = FALSE, na = "")
write.csv(ivy6, paste0(outputLOC, "cdc_ivy6_", today, ".csv"), row.names = FALSE, na = "")
write.csv(ivy7, paste0(outputLOC, "cdc_ivy7_", today, ".csv"), row.names = FALSE, na = "")


#table(seq_list$pangolin_lineage)

#table(ivy5$platename)
