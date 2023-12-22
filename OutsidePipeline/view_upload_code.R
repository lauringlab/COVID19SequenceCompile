################################################################################
#          Creation of VIEW Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 06/23/2023                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################

starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/VIEW_UPLOADS/")

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

seq_list <- filter(seq_list, received_source == "VIEW")


if (length(unique(seq_list$sample_id)) != nrow(seq_list)){
  stop("Duplicate sample IDs - handle accordingly")
}

################################################################################
### fix date formatting
seq_list <- seq_list %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

today <- gsub("-", "", Sys.Date())

view1 <- seq_list %>% select(sample_id, subject_id, coll_date,                    
                            flag, received_source, SiteName, SampleBarcode,                
                            PlateDate, PlatePlatform, PlateNumber,                 
                            pangolin_lineage, pangolin_probability, pangolin_status,             
                            pangolin_note, nextclade_clade, nextclade_totalMissing,      
                            nextclade_completeness, gisaid_strain, gisaid_epi_isl,  
                            genbank_SequenceID, genbank_Accession, 
                            genbank_SubmissionID,
                            received_date, position,          
                            PlateName, PlatePosition, SampleSourceLocation,        
                            pangoLEARN_version, pango_version,               
                            pangolin_version, pangolin_conflict,     
                            nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,    
                            nextclade_totalNonACGTNs, 
                            data_quality_rule, newest_pangolin_lineage, newest_pangolin_date, sample_id_lauring)

colnames(view1) <- c("sample_id", "subject_id", "coll_date",                    
                    "flag", "received_source", "sitename", "samplebarcode",                
                    "platedate", "plateplatform", "platenumber",                 
                    "pangolin_lineage", "pangolin_probability", "pangolin_status",             
                    "pangolin_note", "nextclade_clade", "nextclade_totalmissing",      
                    "nextclade_completeness", "gisaid_strain", "gisaid_epi_isl",
                    "genbank_sequenceid", "genbank_accession", 
                    "genbank_submissionid",
                    "received_date", "position",          
                    "platename", "plateposition", "samplesourcelocation",        
                    "pangolearn_version", "pango_version",               
                    "pangolin_version", "pangolin_conflict",     
                    "nextclade_qcoverallscore", "nextclade_qcoverallstatus", "nextclade_totalmutations",    
                    "nextclade_totalnonacgtns", 
                    "data_quality_rule", "newest_pangolin_lineage", "newest_pangolin_date", "public_gisaid")


#seq_list <- filter(seq_list, platenumber <= 49)
write.csv(view1, paste0(outputLOC, "view_", today, ".csv"), row.names = FALSE, na = "")

#table(seq_list$pangolin_lineage)
