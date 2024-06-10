################################################################################
#          Creation of RIGHT Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 06/03/2024                             #
#                 Code Edited By: Leigh Papalambros                       #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################
checking_wd <- getwd()

if (grepl("leighbak", checking_wd)){
  starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"
  outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/")
  
}else if (grepl("chbl", checking_wd)){
  starting_path <- "/Users/chbl/Dropbox (University of Michigan)/MED-LauringLab/"
  outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/")
} else {
  
  print("User not recognized.")
  
}

seq_list <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################
# 
# seq_list <- seq_list %>% select(sample_id, subject_id, coll_date,                    
#                                 flag, received_source, SiteName, SampleBarcode,                
#                                 PlateDate, PlatePlatform, PlateNumber,                 
#                                 pangolin_lineage, pangolin_probability, pangolin_status,             
#                                 pangolin_note, nextclade_clade, nextclade_totalMissing,      
#                                 nextclade_completeness, gisaid_strain, gisaid_epi_isl,
#                                 genbank_sequenceid, genbank_accession, 
#                                 genbank_submissionid,
#                                 received_date, position,          
#                                 PlateName, PlatePosition, SampleSourceLocation,        
#                                 pangoLEARN_version, pango_version,               
#                                 pangolin_version, pangolin_conflict,     
#                                 nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,    
#                                 nextclade_totalNonACGTNs, 
#                                 data_quality_rule, newest_pangolin_lineage, newest_pangolin_date, sample_id_lauring)

seq_list <- filter(seq_list, received_source == "RIGHT")


if (length(unique(seq_list$sample_id)) != nrow(seq_list)){
  stop("Duplicate sample IDs - handle accordingly")
}

################################################################################
### fix date formatting
seq_list <- seq_list %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

today <- gsub("-", "", Sys.Date())

right1 <- seq_list %>% select(sample_id, subject_id, coll_date,                    
                            flag, received_source, received_date, SiteName, SampleBarcode,                
                            PlateName, PlateDate, PlatePlatform, PlateNumber, PlatePosition,
                            SampleSourceLocation,
                            pangolin_lineage, pangolin_status, pangolin_note,
                            pangolin_version, pangolin_conflict, nextclade_clade,
                            nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallScore,
                            nextclade_qcOverallStatus, nextclade_totalMutations, nextclade_totalNonACGTNs,
                            genbank_SequenceID, genbank_Accession, genbank_SubmissionID, data_quality_rule,
                            newest_pangolin_lineage, newest_pangolin_date)
                            

colnames(right1) <- c("sample_id", "subject_id", "coll_date",                    
                    "flag", "received_source", "received_date", "sitename", "samplebarcode",                
                    "platename", "platedate", "plateplatform", "platenumber", "plateposition", 
                    "samplesourcelocation",
                    "pangolin_lineage", "pangolin_status","pangolin_note",
                    "pangolin_version", "pangolin_conflict", "nextclade_clade",
                    "nextclade_totalmissing", "nextclade_completeness", "nextclade_qcoverallscore",
                    "nextclade_qcoverallstatus", "nextclade_totalsubstitutions", "nextclade_totalnonacgtns", 
                    "genbank_sequenceid", "genbank_accession", "genbank_submissionid", "data_quality_rule",
                    "newest_pangolin_lineage", "newest_pangolin_date")
                    



#seq_list <- filter(seq_list, platenumber <= 49)
write.csv(right1, paste0(outputLOC, "right_", today, ".csv"), row.names = FALSE, na = "")

#table(seq_list$pangolin_lineage)
