################################################################################
#       Creation of CDC IVY Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 05/18/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
outputLOC <- paste0(starting_path, "SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/")

seq_list <- read.csv(paste0(starting_path, "SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

seq_list <- seq_list %>% select(sample_id, subject_id, coll_date,                    
                                flag, received_source, SiteName, SampleBarcode,                
                                PlateDate, PlatePlatform, PlateNumber,                 
                                pangolin_lineage, pangolin_probability, pangolin_status,             
                                pangolin_note, nextclade_clade, nextclade_totalMissing,      
                                nextclade_completeness, gisaid_strain, gisaid_epi_isl,              
                                received_date, position, subject_id_length,           
                                PlateName, PlatePosition, SampleSourceLocation,        
                                pangoLEARN_version, pangolin_conflict, pango_version,               
                                pangolin_version, pangolin_runDate, PlateToPangolin_days,        
                                nextclade_qcOverallScore, nextclade_qcOverallStatus, nextclade_totalMutations,    
                                nextclade_totalNonACGTNs, nextclade_runDate, PlateToNextclade_days)

seq_list <- filter(seq_list, received_source == "CDCIVY")

# only keep complete rows?

################################################################################
# After the first upload, we'll need to keep track of what has already been uploaded
# so we'll read in the full list, then read in the previous upload list, and only keep 
# rows that are not in the previous upload list(s) to write out and upload the next time

colnames(seq_list) <- tolower(colnames(seq_list))

# add leading zero to month
if (length(month(Sys.Date()))){
  m <- paste0("0", month(Sys.Date()))
} else {
  m <- month(Sys.Date())
}
# add leading zero to day
if (length(day(Sys.Date()))){
  d <- paste0("0", day(Sys.Date()))
} else {
  d <- day(Sys.Date())
}

today <- paste0(year(Sys.Date()), m, d)

write.csv(seq_list, paste0(outputLOC, "cdc_ivy_", today, ".csv"), row.names = FALSE, na = "")
