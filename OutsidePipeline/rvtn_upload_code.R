################################################################################
#          Creation of RVTN Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 03/07/2022                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################

starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/RVTN_UPLOADS/")

seq_list <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

################################################################################

seq_list <- seq_list %>% select(sample_id, subject_id, coll_date,                    
                                flag, received_source, SiteName, SampleBarcode,                
                                PlateDate, PlatePlatform, PlateNumber,                 
                                pangolin_lineage, pangolin_probability, pangolin_status,             
                                pangolin_note, nextclade_clade, nextclade_totalMissing,      
                                nextclade_completeness, gisaid_strain, gisaid_epi_isl,              
                                received_date, position,          
                                PlateName, PlatePosition, SampleSourceLocation,        
                                pangoLEARN_version, pangolin_conflict, pango_version,               
                                pangolin_version, pangolin_runDate, PlateToPangolin_days,        
                                nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,    
                                nextclade_totalNonACGTNs, nextclade_runDate, PlateToNextclade_days)

seq_list <- filter(seq_list, received_source == "RVTN")


if (length(unique(seq_list$sample_id)) != nrow(seq_list)){
  stop("Duplicate sample IDs - handle accordingly")
}

################################################################################
### fix date formatting
seq_list <- seq_list %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))


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

rvtn <- seq_list %>% select(sample_id, subject_id, coll_date,                    
                                flag, received_source, sitename, samplebarcode,                
                                platedate, plateplatform, platenumber, 
                                received_date, position, platename, plateposition, samplesourcelocation,
                                gisaid_strain, gisaid_epi_isl, pangolearn_version,
                                pango_version, pangolin_version, pangolin_lineage, pangolin_status,             
                                pangolin_note, nextclade_clade, nextclade_totalmissing,      
                                nextclade_completeness, nextclade_qcoverallscore, nextclade_qcoverallstatus, 
                                nextclade_totalmutations, pangolin_probability, nextclade_totalnonacgtns)


#seq_list <- filter(seq_list, platenumber <= 49)
write.csv(rvtn, paste0(outputLOC, "rvtn_", today, ".csv"), row.names = FALSE, na = "")

#table(seq_list$pangolin_lineage)