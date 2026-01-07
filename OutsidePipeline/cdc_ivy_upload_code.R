################################################################################
#       Creation of CDC IVY Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 07/20/2023                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)
library(janitor)
library(dplyr)
library(stringr)

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

samp_dup <- seq_list %>% group_by(sample_id) %>% summarize(count = n ()) %>% filter(count > 1)
print(samp_dup)
# these samples dupilcated because they have different "newest_nextclade_clades"
dup_samples <- filter(seq_list, sample_id %in% c("G43Q32V2", "W13J94M7", "W13J95C6", "W13J95O2"))


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
ivy7 <- filter(seq_list, substr(subject_id, 1, 2) == 25)


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

ivy6$newest_nextclade_date <-ymd(ivy6$newest_nextclade_date)


ivy6 <- ivy6 %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                        received_date, sitename, samplebarcode, platename, platedate, 
                        plateplatform, platenumber, plateposition, samplesourcelocation, 
                        pangolin_lineage, pangolin_status, pangolin_note, 
                        pangolin_version, pangolin_conflict, nextclade_clade,
                        nextclade_totalmissing, nextclade_completeness, nextclade_qcoverallscore,
                        nextclade_qcoverallstatus, nextclade_totalmutations, 
                        nextclade_totalnonacgtns, genbank_sequenceid, genbank_accession, 
                        genbank_submissionid, data_quality_rule, newest_pangolin_lineage, 
                        newest_pangolin_date, sf456l_present, newest_nextclade_clade, newest_nextclade_date)

colnames(ivy6) <- c("sample_id", "subject_id", "coll_date", "flag", "received_source", 
                    "received_date", "sitename", "samplebarcode", "platename", "platedate", 
                    "plateplatform", "platenumber", "plateposition", "samplesourcelocation", 
                    "pangolin_lineage", "pangolin_status", "pangolin_note", 
                    "pangolin_version", "pangolin_conflict", "nextclade_clade",
                    "nextclade_totalmissing", "nextclade_completeness", "nextclade_qcoverallscore",
                    "nextclade_qcoverallstatus", "nextclade_totalsubstitutions", 
                    "nextclade_totalnonacgtns", "genbank_sequenceid", "genbank_accession", 
                    "genbank_submissionid", "data_quality_rule", "newest_pangolin_lineage", 
                    "newest_pangolin_date", "nextclade_sf456l_present", "newest_nextclade_clade", "newest_nextclade_date")

################################################################################
## for IVY7 we are including an over write of the nextclade clades to the newest assignments
## (ie like we do for the updating of pangolin) but they are not new columns 


ivy7$newest_nextclade_date <-ymd(ivy7$newest_nextclade_date)
ivy7$newest_pangolin_date <-ymd(ivy7$newest_pangolin_date)
ivy7$sample_id <-trimws(ivy7$sample_id)
ivy7$subject_id <-trimws(ivy7$subject_id)


ivy7 <- ivy7 %>% select(sample_id, subject_id, coll_date, flag, received_source, 
                        received_date, sitename, samplebarcode, platename, platedate, 
                        plateplatform, platenumber, plateposition, samplesourcelocation, 
                        pangolin_lineage, pangolin_status, pangolin_note, 
                        pangolin_version, pangolin_conflict, nextclade_clade,
                        nextclade_totalmissing, nextclade_completeness, nextclade_qcoverallscore,
                        nextclade_qcoverallstatus, nextclade_totalmutations, 
                        nextclade_totalnonacgtns, genbank_sequenceid, genbank_accession, 
                        genbank_submissionid, data_quality_rule, newest_pangolin_lineage, 
                        newest_pangolin_date, sf456l_present, newest_nextclade_clade, newest_nextclade_date)

colnames(ivy7) <- c("sample_id", "subject_id", "coll_date", "flag", "received_source", 
                    "received_date", "sitename", "samplebarcode", "platename", "platedate", 
                    "plateplatform", "platenumber", "plateposition", "samplesourcelocation", 
                    "pangolin_lineage", "pangolin_status", "pangolin_note", 
                    "pangolin_version", "pangolin_conflict", "nextclade_clade",
                    "nextclade_totalmissing", "nextclade_completeness", "nextclade_qcoverallscore",
                    "nextclade_qcoverallstatus", "nextclade_totalsubstitutions", 
                    "nextclade_totalnonacgtns", "genbank_sequenceid", "genbank_accession", 
                    "genbank_submissionid", "data_quality_rule", "newest_pangolin_lineage", 
                    "newest_pangolin_date", "nextclade_sf456l_present", "newest_nextclade_clade", "newest_nextclade_date")

#seq_list <- filter(seq_list, platenumber <= 49)
#write.csv(ivy3, paste0(outputLOC, "cdc_ivy3_", today, ".csv"), row.names = FALSE, na = "")
#write.csv(ivy4, paste0(outputLOC, "cdc_ivy4_", today, ".csv"), row.names = FALSE, na = "")
#write.csv(ivy5, paste0(outputLOC, "cdc_ivy5_", today, ".csv"), row.names = FALSE, na = "")
#write.csv(ivy6, paste0(outputLOC, "cdc_ivy6_", today, ".csv"), row.names = FALSE, na = "")
write.csv(ivy7, paste0(outputLOC, "cdc_ivy7_", today, ".csv"), row.names = FALSE, na = "")


#table(seq_list$pangolin_lineage)

#table



################################################################################

## Custom processing of IVY7 KP.2 spike mutations ##

## Wuhan reference section
#ivy7_nc_p1 <- read_tsv(paste0(nc_fp, "ivy7_1_nextclade.tsv"))
#dels <- ivy7_nc_p1 %>% select(aaDeletions)

## Read in the nextclade.tsv files for all IVY7 samples
nc_fp <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/External_Projects_DataRequests/SARSCOV2/kp2/"

file_list <- list.files(pattern = "*_nextclade.tsv", path = nc_fp)

nc_storage <- data.frame()

for (each_page in file_list){
  nc1 <- read.table(paste0(nc_fp, "/", each_page), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
  
  #print(colnames(nc1))
  if("totalMutations" %in% colnames(nc1)){
    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalMutations, totalNonACGTNs, aaSubstitutions, aaDeletions)
  } else {
    ### nextclade update changed column totalMutations to totalSubstitutions (near 6/18/2021)
    nc1 <- nc1 %>% select(seqName, clade, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs, aaSubstitutions, aaDeletions)
    colnames(nc1) <- c("seqName", "clade", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs", "aaSubstitutions", "aaDeletions")
  }
  #nc1$nextclade_runDate <- date_from_file_FIRST(each_page)
  
  ## the KP.2 spike mutations substitions only of interest (55 total) + 5 deletions
  nc1 <- nc1 %>% mutate(ST19I_present = case_when(grepl("S:T19I", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SR21T_present = case_when(grepl("S:R21T", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SL24S_present = case_when(grepl("S:L24S", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SS50L_present = case_when(grepl("S:S50L", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SV127F_present = case_when(grepl("S:V127F", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SG142D_present = case_when(grepl("S:G142D", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SF157S_present = case_when(grepl("S:F157S", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SR158G_present = case_when(grepl("S:R158G", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN211I_present = case_when(grepl("S:N211I", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SV213G_present = case_when(grepl("S:V213G", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SL216F_present = case_when(grepl("S:L216F", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SH245N_present = case_when(grepl("S:H245N", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SA264D_present = case_when(grepl("S:A264D", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SI332V_present = case_when(grepl("S:I332V", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SG339H_present = case_when(grepl("S:G339H", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SR346T_present = case_when(grepl("S:R346T", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SK356T_present = case_when(grepl("S:K356T", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SS371F_present = case_when(grepl("S:S371F", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SS373P_present = case_when(grepl("S:S373P", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SS375F_present = case_when(grepl("S:S375F", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(ST376A_present = case_when(grepl("S:T376A", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SR403K_present = case_when(grepl("S:R403K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SD405N_present = case_when(grepl("S:D405N", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SR408S_present = case_when(grepl("S:R408S", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SK417N_present = case_when(grepl("S:K417N", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN440K_present = case_when(grepl("S:N440K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SV445H_present = case_when(grepl("S:V445H", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SG446S_present = case_when(grepl("S:G446S", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN450D_present = case_when(grepl("S:N450D", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SL452W_present = case_when(grepl("S:L452W", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SL455S_present = case_when(grepl("S:L455S", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SF456L_present = case_when(grepl("S:F456L", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN460K_present = case_when(grepl("S:N460K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SS477N_present = case_when(grepl("S:S477N", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(ST478K_present = case_when(grepl("S:T478K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN481K_present = case_when(grepl("S:N481K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SE484K_present = case_when(grepl("S:E484K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SF486P_present = case_when(grepl("S:F486P", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SQ498R_present = case_when(grepl("S:Q498R", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN501Y_present = case_when(grepl("S:N501Y", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SY505H_present = case_when(grepl("S:Y505H", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SE554K_present = case_when(grepl("S:E554K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SA570V_present = case_when(grepl("S:A570V", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SD614G_present = case_when(grepl("S:D614G", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SP621S_present = case_when(grepl("S:P621S", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SH655Y_present = case_when(grepl("S:H655Y", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN679K_present = case_when(grepl("S:N679K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SP681R_present = case_when(grepl("S:P681R", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN764K_present = case_when(grepl("S:N764K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SD796Y_present = case_when(grepl("S:D796Y", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SS939F_present = case_when(grepl("S:S939F", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SQ954H_present = case_when(grepl("S:Q954H", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SN969K_present = case_when(grepl("S:N969K", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SV1104L_present = case_when(grepl("S:V1104L", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(SP1143L_present = case_when(grepl("S:P1143L", aaSubstitutions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(Sdel25_present = case_when(grepl("S:P25", aaDeletions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(Sdel69_present = case_when(grepl("S:H69", aaDeletions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(Sdel70_present = case_when(grepl("S:V70", aaDeletions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(Sdel144_present = case_when(grepl("S:Y144", aaDeletions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(Sdel212_present = case_when(grepl("S:N211", aaDeletions) ~ 1, T ~ 0))
  nc1 <- nc1 %>% mutate(Sdel483_present = case_when(grepl("S:V483", aaDeletions) ~ 1, T ~ 0))
  
  
  nc_storage <- rbind(nc_storage, nc1)
}

## just looking at the aa substituions column

only_aa_subs <- nc1 %>% select(aaSubstitutions)


## mathcing just the IVY7 samples from the nextclade file
## select only the IVY7 sample_ids from the nextclade file

nc_storage <- nc_storage %>%
  mutate(total_aa_subs = rowSums(select(., 10:64)))

nc_storage <- nc_storage %>%
  mutate(total_del_subs = rowSums(select(., 65:70)))

nc_storage <- nc_storage %>%
  mutate(total_difference_aa_subs_wuhan = 55 - total_aa_subs)

#nc_storage$SaaSubstitutions <- nc_storage %>% select(grepl("S:"), aaSubstitutions)
#pattern = "^,S:"
#nc_storage <- nc_storage %>%
#  mutate(SaaSubstitutions = ifelse(grepl(pattern, aaSubstitutions), aaSubstitutions, NA))

#select only the new nextclade columns

nc_storage_small <- nc_storage %>% select(seqName, aaSubstitutions, aaDeletions, 10:73)

i7_w_muts <- merge(nc_storage_small, ivy7, by.x = "seqName", by.y = "sample_id", all.y = TRUE)



#####################################################################################
## the BA.2.86 nextclade reference version section

#one of the nextclade files just to look at it raw
ivy7_nc_ba2_p1 <- read_tsv(paste0(nc_fp, "ivy7_1_ba2nextclade.tsv"))
#colnames(ivy7_nc_ba2_p1)

file_list_2 <- list.files(pattern = "*_ba2nextclade.tsv", path = nc_fp)

nc_storage_2 <- data.frame()

for (each_page in file_list_2){
  nc2 <- read.table(paste0(nc_fp, "/", each_page), header = TRUE, colClasses = "character", stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
  
  #print(colnames(nc1))
  if("totalMutations" %in% colnames(nc2)){
    nc2 <- nc2 %>% select(seqName, clade, clade_display, totalMissing, qc.overallScore, qc.overallStatus, totalMutations, totalNonACGTNs, aaSubstitutions, aaDeletions)
  } else {
    ### nextclade update changed column totalMutations to totalSubstitutions (near 6/18/2021)
    nc2 <- nc2 %>% select(seqName, clade, clade_display, totalMissing, qc.overallScore, qc.overallStatus, totalSubstitutions, totalNonACGTNs, aaSubstitutions, aaDeletions, 56, 57)
    colnames(nc2) <- c("seqName", "clade", "clade_display", "totalMissing", "qc.overallScore", "qc.overallStatus", "totalMutations", "totalNonACGTNs", "aaSubstitutions", "aaDeletions", "aaSubstitutions_JN1", "aaDeletions_JN1")
  }

  nc2 <- nc2 %>% mutate(SL455S_present = case_when(grepl("S:L455S", aaSubstitutions) ~ 1, T ~ 0))
  nc2 <- nc2 %>% mutate(SR346T_present = case_when(grepl("S:R346T", aaSubstitutions_JN1) ~ 1, T ~ 0))
  nc2 <- nc2 %>% mutate(SF456L_present = case_when(grepl("S:F456L",aaSubstitutions_JN1) ~ 1, T ~ 0))
  nc2 <- nc2 %>% mutate(SV1104L_present = case_when(grepl("S:V1104L", aaSubstitutions_JN1) ~ 1, T ~ 0))
  nc2 <- nc2 %>% mutate(S31_del_present = case_when(grepl("S:S31-", aaDeletions_JN1) ~ 1, T ~ 0))
  
  
  nc_storage_2 <- rbind(nc_storage_2, nc2)
}


#print(file_list_2)
## pulling out just the spike information
nc_storage_2$just_spike_aasubs <- str_extract_all(nc_storage_2$aaSubstitutions_JN1, "S:[A-Z0-9]{4,6}")
nc_storage_2$just_spike_aadels <- str_extract_all(nc_storage_2$aaDeletions_JN1, "S:[A-Z0-9-]{4,6}")
#print(unlist(nc2$spike_muts))


nc_storage_2$just_spike_aasubs <- sapply(nc_storage_2$just_spike_aasubs, paste, collapse = ", ")
nc_storage_2$just_spike_aadels <- sapply(nc_storage_2$just_spike_aadels, paste, collapse = ", ")


## list of KP.2 aaSubstitutions and deletions
kp2_muts <- "S:R346T|S:F456L|S:V1104L"

nc_storage_2$non_kp2_spike_aasubs <- str_remove_all(nc_storage_2$just_spike_aasubs, kp2_muts)

nc_storage_2$non_kp2_spike_aasubs <- sapply(nc_storage_2$non_kp2_spike_aasubs, paste, collapse = c(",", ", ", ""))

nc_storage_2$non_kp2_spike_aasubs <- gsub(",\\s+", " ", nc_storage_2$non_kp2_spike_aasubs)
nc_storage_2$non_kp2_spike_aasubs <- trimws(nc_storage_2$non_kp2_spike_aasubs)

kp2_dels <- "S:S31-"
nc_storage_2$non_kp2_spike_aadels <- str_remove_all(nc_storage_2$just_spike_aadels, kp2_dels)
nc_storage_2$non_kp2_spike_aadels <- gsub(",\\s+", " ", nc_storage_2$non_kp2_spike_aadels)
nc_storage_2$non_kp2_spike_aadels <- trimws(nc_storage_2$non_kp2_spike_aadels)


##################
## looking at KP.2 aa subs counts
# looking at other versions of S:R346 that isn't the KP.2 specific version

#nc2$non_kp2_spike <- str_trim(nc2$non_kp2_spike)

#nc2$spike_muts <- as.character(nc2$spike_muts)

#nc2$non_kp2_spike <- nc2$non_kp2_spike[nc2$non_kp2_spike != ""]

#str(nc2$spike_muts)


## a function to count pattern occurrences of S:**** and S:*****

#count_patterns <- function(text) {
#  count_s_char <- str_count(text, "S:")
  #count_five_char <- str_count(text, "S:[A-Z0-9]{5}")
  #count_six_char <- str_count(text, "S:[A-Z0-9]{6}")
  #four_count <- count_four_char
  #five_count <- count_five_char
  #six_count <- count_six_char
#  total_count <- count_s_char
#  return(total_count)
#}

count_patterns <- function(text) {
  count_s_char <- str_count(text, "S:")
  total_count <- count_s_char
  return(total_count)
}


nc_storage_2 <- nc_storage_2 %>%
  mutate(num_spike_aasubs_nonkp2 = sapply(non_kp2_spike_aasubs, count_patterns))

nc_storage_2 <- nc_storage_2 %>%
  mutate(num_spike_aadels_nonkp2 = sapply(non_kp2_spike_aadels, count_patterns))


## this section is creating the columns for S1, RBD, and NTD mutation counts

## S1 aa range is 14-685
number_range_pattern_s1 <- "\\bS:[A-Z](1[4-9]|[2-9][0-9]|[1-5][0-9]{2}|6[0-7][0-9]|68[0-5])[A-Z]\\b"

test_pull <- nc_storage_2 %>%
  mutate(s1_aa_subs_spike_only_nonkp2 = sapply(
    str_extract_all(nc_storage_2$non_kp2_spike_aasubs, number_range_pattern_s1),
    paste, collapse = " "
  ))

test_pull <- test_pull %>%
  mutate(num_spike_aasubs_nonkp2_S1 = sapply(s1_aa_subs_spike_only_nonkp2, count_patterns))

## RBD range 319-541
number_range_pattern_RBD <- "\\bS:[A-Z](3[1-9][0-9]|4[0-9][0-9]|5[0-3][0-9]|54[0-1])[A-Z]\\b"

test_pull_2 <- test_pull %>%
  mutate(rbd_aa_subs_spike_only_nonkp2 = sapply(
    str_extract_all(test_pull$non_kp2_spike_aasubs, number_range_pattern_RBD),
    paste, collapse = " "
  ))

test_pull_2 <- test_pull_2 %>%
  mutate(num_spike_aasubs_nonkp2_RBD = sapply(rbd_aa_subs_spike_only_nonkp2, count_patterns))

## NTD range 14-305
number_range_pattern_NTD <- "\\bS:[A-Z](1[4-9]|[2-9][0-9]|1[0-9][0-9]|2[0-9][0-9]|30[0-5])[A-Z]\\b"

test_pull_3 <- test_pull_2 %>%
  mutate(ntd_aa_subs_spike_only_nonkp2 = sapply(
    str_extract_all(test_pull_2$non_kp2_spike_aasubs, number_range_pattern_NTD),
    paste, collapse = " "
  ))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_NTD = sapply(ntd_aa_subs_spike_only_nonkp2, count_patterns))
#nc_storage_2 <- nc_storage_2 %>%
#  mutate(total_aa_subs = rowSums(select(., 12:15)))

#nc_storage_2 <- nc_storage_2 %>%
 # mutate(total_differences_ba2 = 4 - total_aa_subs)


## need to add a piece of code looking for S:-31 in the aasubs column (ancestral)
# if there is a "0" in the SR346T_present, SF456L_present, and SV1104L_present columns then we need to add
# a 1 to the total count of non kp spike aa subs and the corresponding domain regions

# accounting for the ancestral changes of SR346T_present in the total AAsubs and the RBD and S1 domains
test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2 = if_else(SR346T_present == 0, num_spike_aasubs_nonkp2 + 1, num_spike_aasubs_nonkp2))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_RBD = if_else(SR346T_present == 0, num_spike_aasubs_nonkp2_RBD + 1, num_spike_aasubs_nonkp2_RBD))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_S1 = if_else(SR346T_present == 0, num_spike_aasubs_nonkp2_S1 + 1, num_spike_aasubs_nonkp2_S1))

# accounting for the ancestral changes of SF456L_present in the total AAsubs and the RBD and S1 domains

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2 = if_else(SF456L_present == 0, num_spike_aasubs_nonkp2 + 1, num_spike_aasubs_nonkp2))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_RBD = if_else(SF456L_present == 0, num_spike_aasubs_nonkp2_RBD + 1, num_spike_aasubs_nonkp2_RBD))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_S1 = if_else(SF456L_present == 0, num_spike_aasubs_nonkp2_S1 + 1, num_spike_aasubs_nonkp2_S1))

# accounting for the ancestral changes of SV1104L_present in the total AA subs

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2 = if_else(SV1104L_present == 0, num_spike_aasubs_nonkp2 + 1, num_spike_aasubs_nonkp2))


# accounting for the ancestral deletion of S31_del_present in the total AAsubs and the S1 and NTD domains

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2 = if_else(S31_del_present == 1, num_spike_aasubs_nonkp2 + 1, num_spike_aasubs_nonkp2))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_NTD = if_else(S31_del_present == 1, num_spike_aasubs_nonkp2_NTD + 1, num_spike_aasubs_nonkp2_NTD))

test_pull_3 <- test_pull_3 %>%
  mutate(num_spike_aasubs_nonkp2_S1 = if_else(S31_del_present == 1, num_spike_aasubs_nonkp2_S1 + 1, num_spike_aasubs_nonkp2_S1))


nc_storage_2_small <- test_pull_3 %>% select(seqName, clade_display, aaSubstitutions, aaDeletions, aaSubstitutions_JN1, aaDeletions_JN1, SL455S_present,
                                              SR346T_present, SF456L_present, SV1104L_present, S31_del_present, just_spike_aasubs, just_spike_aadels, non_kp2_spike_aasubs,
                                              num_spike_aasubs_nonkp2, non_kp2_spike_aadels, num_spike_aadels_nonkp2, s1_aa_subs_spike_only_nonkp2,
                                             num_spike_aasubs_nonkp2_S1, rbd_aa_subs_spike_only_nonkp2, num_spike_aasubs_nonkp2_RBD, 
                                             ntd_aa_subs_spike_only_nonkp2, num_spike_aasubs_nonkp2_NTD )
#str(nc_storage_2_small)

i7_ba_muts <- merge(nc_storage_2_small, ivy7, by.x = "seqName", by.y = "sample_id", all.y = TRUE)

#i7_ba_muts <- apply(i7_ba_muts, as.character, FUN = paste)
#nc_storage_2$non_kp2_spike_aasubs <- vapply(nc_storage_2$non_kp2_spike_aasubs, paste, collapse = ", ", character(1L))

#################################################################################
##Cleaning up the data and formatting it section

## trying to break the list of spike mutations that are left into new columns
#df_expanded <- nc_storage_2 %>%
#  unnest_wider(non_kp2_spike_aasubs, names_sep = ",")

#duplicating the list of non spike aaSubs as a back up/double check
nc_storage_2_small$non_kp2_spike_aasubs_2_2 <- nc_storage_2_small$non_kp2_spike_aasubs_2

df_separated <- nc_storage_2_small %>%
  separate(non_kp2_spike_aasubs_2, into = c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7", "Col8", "Col9", "Col10", "Col11", "Col12"), sep = ", ", fill = "right")

## putting just the values that have spike aasubs into a new column 
df_concat <- df_separated %>%
  rowwise() %>%
  mutate(non_kp_only_spike_aasubs = paste(na.omit(c_across(Col1:Col12)[nzchar(c_across(Col1:Col12))]), collapse = " "))

## formatting num_spike_aadels


#################################################################################




##############################

## just selecting the non redcap information for a write file

wf_i7_ba_muts <- i7_ba_muts %>% select(seqName, subject_id, clade_display, aaSubstitutions, aaDeletions, aaSubstitutions_JN1, aaDeletions_JN1, SL455S_present,
                                       SR346T_present, SF456L_present, SV1104L_present, S31_del_present, just_spike_aasubs, just_spike_aadels, non_kp2_spike_aasubs,
                                       num_spike_aasubs_nonkp2, non_kp2_spike_aadels, num_spike_aadels_nonkp2, s1_aa_subs_spike_only_nonkp2,
                                       num_spike_aasubs_nonkp2_S1, rbd_aa_subs_spike_only_nonkp2, num_spike_aasubs_nonkp2_RBD, 
                                       ntd_aa_subs_spike_only_nonkp2, num_spike_aasubs_nonkp2_NTD)


df_Place2 = data.frame(lapply(wf_i7_ba_muts, as.character), stringsAsFactors=FALSE)

#wf_i7_wu_muts <- i7_w_muts %>% select(seqName, aaSubstitutions, aaDeletions, 4:67)

write.csv(wf_i7_ba_muts, paste0(nc_fp, "IVY7_BA286_spike_6.csv"))


write.csv(wf_i7_wu_muts, paste0(nc_fp, "IVY7_WUHAN_muts.csv"))


#last_upload <- read.csv("~/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/cdc_ivy7_20250530.csv")

kevins_missing_samples <- read.csv("~/University of Michigan Dropbox/MED-LauringLab/External_Projects_DataRequests/SARSCOV2/kp2/2025-06-24_90-missing-samples.csv")

# no matchs to our ivy7 data under sample_id but under subject_id we have matches
our_data_kevins_samples <- merge(kevins_missing_samples, ivy7, by.x = "aliquot_id", by.y = "sample_id", all.x = TRUE)

our_data_kevins_samples_2 <- merge(kevins_missing_samples, ivy7, by.x = "study_id", by.y = "subject_id", all.x = TRUE)

# counting the instances of subject_id to make sure all of them match the correct sample_id

sub_id_totals <- wf_i7_ba_muts %>% group_by(subject_id) %>% summarize(count = n ())

dups_sub_id <- filter(sub_id_totals, count > 1)

#list of the duplicate subject_ids and the corresponding sample/aliquot_id

dups_sub_id_list <- merge(dups_sub_id, wf_i7_ba_muts, by = "subject_id", all.x = TRUE)

write.csv(dups_sub_id_list, paste0(nc_fp, "dup_subject_id_BA286.csv"))

