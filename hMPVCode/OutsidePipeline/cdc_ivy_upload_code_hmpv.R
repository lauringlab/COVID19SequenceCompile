### Code to generate RedCap upload file for hMPV - IVY


library(tidyverse)
library(lubridate)


# read in full compiled hMPV file

hMPV <- read.csv("/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")

ivy_hmpv <- hMPV %>% filter(grepl("IVY", received_source))



#sapply(rsv_ab_out5, class)
#ivy_hmpv$genbank_SubmissionID <- ""
#ivy_hmpv$genbank_Accession <- ""
#ivy_hmpv$genbank_SequenceID <- ""

hmpv_out <- ivy_hmpv %>% select(subject_id, sample_id, coll_date, flag, 
                                  received_source, received_date, SiteName, SampleBarcode,
                                  PlateName, PlateDate, PlatePlatform, PlateNumber, PlatePosition, SampleSourceLocation, 
                                  nextclade_clade, nextclade_alternate_clade, nextclade_alternate_clade, nextclade_totalMissing, nextclade_completeness, nextclade_qcOverallScore,
                                nextclade_qcOverallStatus, nextclade_totalMutations, nextclade_totalNonACGTNs, 
                                genbank_SequenceID, genbank_Accession, genbank_SubmissionID, data_quality_rule)




colnames(hmpv_out) <- c("subject_id", "sample_id", "coll_date_hmpv", "flag_hmpv", 
                           "received_source_hmpv", "received_date_hmpv", "sitename_hmpv", "samplebarcode_hmpv",
                           "plate_name_hmpv", "plate_date_hmpv", "plate_platform_hmpv", "plate_number_hmpv", "plate_position_hmpv", "sample_source_location_hmpv",
                           "nextclade_clade_hmpv", "nextclade_alternate_clade_hmpv", "nextclade_totalmissing_hmpv", "nextclade_completeness_hmpv", "nextclade_qcoverallscore_hmpv",
                        "nextclade_qcoverallstatus_hmpv", "nextclade_totalsubstitutions_hmpv", "nextclade_totalnonacgtns_hmpv", 
                        "genbank_sequenceid_hmpv", "genbank_accession_hmpv", "genbank_submissionid_hmpv", "data_quality_rule_hmpv")


f_out <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/SEQUENCING/hMPV/4_SequenceSampleMetadata/FinalSummary/IVY_uploads/"



write.csv(hmpv_out, paste0(f_out, "ivy_hmpv_upload_", gsub("-", "", Sys.Date()), ".csv"), row.names = FALSE, na = "")




