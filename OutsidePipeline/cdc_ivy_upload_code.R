################################################################################
#       Creation of CDC IVY Upload Dataset for COVID-19 Genetic Sampling       #
#                         Last Updated: 06/18/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/CDC_IVY_UPLOADS/")

seq_list <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")
seq_list_o <- read.csv(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character")

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

seq_list_o <- seq_list_o %>% select(sample_id, subject_id, coll_date,                    
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


seq_list <- filter(seq_list, received_source == "CDCIVY")

if (length(unique(seq_list$sample_id)) != nrow(seq_list)){
  stop("Duplicate sample IDs - handle accordingly")
}

## remove study withdraws
seq_list <- filter(seq_list, flag != "Withdrawn from study")

################################################################################
### fix date formatting
seq_list <- seq_list %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

seq_list_o <- seq_list_o %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                                      grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                                      T ~ NA_character_))

################################################################################
# create portion for qPCR matching

qpcr <- seq_list %>% select(sample_id, subject_id, coll_date, PlateName, PlatePosition)

qpcr_full <- filter(seq_list_o, PlateName %in% unique(qpcr$PlateName)) %>% select(sample_id, subject_id, coll_date, PlateName, PlatePosition, received_source, flag)

pmc <- read.csv("C:/Users/juliegil/Documents/UofM_Work/Lauring_Lab/plate_map_cross.csv")
qpcr_full <- merge(qpcr_full, pmc, by.x = c("PlatePosition"), by.y = c("Slot"), all.x = TRUE)
qpcr_full$wellgrid <- paste0(qpcr_full$Letter, qpcr_full$Number) 
qpcr_full <- qpcr_full %>% arrange(PlateName, Letter, Number)

write.csv(qpcr_full, "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/8_QPCR_IVY/IVY_Locations/ivy_mapping.csv", row.names = FALSE, na = "")

# id_count <- seq_list %>% group_by(sample_id) %>% summarize(count = length(subject_id))
# seq_list <- merge(seq_list, id_count, by = c("sample_id"), all.x = TRUE)
# 
# dupes <- filter(seq_list, count > 1)
### if we need to manually put in gisaid info

# only keep complete rows?

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

#seq_list <- filter(seq_list, platenumber <= 49)
write.csv(seq_list, paste0(outputLOC, "cdc_ivy_", today, ".csv"), row.names = FALSE, na = "")

#table(seq_list$pangolin_lineage)
