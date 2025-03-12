################################################################################
#                    GISAID File Upload Format Creation                        #
#                         Last Updated: 06/02/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(openxlsx)

################################################################################
# just need some of these functions

checking_wd <- getwd()
if (grepl("juliegil", checking_wd)){
  
  code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
  #code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
  starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  #code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
  #starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  
} else if (grepl("leighbaker", checking_wd)){
  code_path <- "/Users/leighbaker/Documents/Lauring_Lab/COVID19SequenceCompile/"
  starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"

} else if (grepl("leighbak", checking_wd)){
  
  starting_path <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/"
  code_path <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/"
  batch_path <- "/Users/leighbak/Documents/Lauring_Lab/AlertCode"
    
} else {
  
  print("User not recognized.")
  
}


source(paste0(code_path, "pipeline_functions.R"))

################################################################################
### fill in some info

#fill in the plate name below if running this code seperate and not after "full_run_code.R"
#plate_name <- "20250226_RSVB_Nanopore_Run_13"

plate_datef <- strsplit(plate_name, "_")[[1]][1] # plate date in YYYYMMDD format
runtech <- strsplit(plate_name, "_")[[1]][3] # nanopore or illumina, will match "PlatePlatform" options
runnum <- strsplit(plate_name, "_")[[1]][5] # number, will match "PlateNumber" options

################################################################################

# create gisaid directory
dir.create(paste0(starting_path, "SEQUENCING/RSV_B/6_GenBank_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum))

################################################################################

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
source(paste0(code_path, "OutsidePipeline/checking_compiled_files.R"))

# set output path for gisaid upload file
# will need to add appropriate folder name at the end of this path
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_B/6_GenBank_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum, "/")

################################################################################

# read in full compiled pile
finalfileLOC <- paste0(starting_path, "SEQUENCING/RSV_B/4_SequenceSampleMetadata/FinalSummary")
final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")

final_file <- filter(final_file, !grepl("Missing Date in Manifest", flag))

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_completeness) >= 80)

# select run of choice
ff <- filter(ff, PlatePlatform == runtech & PlateNumber == runnum)

table(ff$received_source, useNA = "always")
#ff <- filter(ff, received_source != "BSL3")
################################################################################
# set up alert to duplicate items

if (any(ff$sample_per_subject > 1)){
  print("STOP: Examine this set of GENBANK submissions.")
  stop("There are samples from subject_ids that we've sequenced previously.")
}



#samples_previous <- filter(ff, sample_per_subject > 1) %>% select(subject_id, sample_id, coll_date)
#original_full <- filter(final_file, subject_id %in% unique(samples_previous$subject_id))
# if they are IVYIC duplicates, let them all through regardless
### check if the samples are > 90 days apart from one another - then you can let 
### them through.

#writes the original_full dataframe to a .csv file
#the path will need to be changed to where you want the file to be written to
#write_csv(df, 'path/to/df/export.csv')
#write_csv(original_full, '~/Documents/Leigh_coding_projects/NP168.csv')

### uncomment this portion to remove those samples
### to remove these: 
#ff <- filter(ff, sample_per_subject == 1 | subject_id %in% c("017838404", "030452124", 
#                                                              "032982827", "045586418", 
#                                                              "100422404"))
#ff <- filter(ff, sample_per_subject == 1)
#ff <- filter(ff, subject_id != "101074339")
#ff <- filter(ff, subject_id != "045447388" & subject_id != "017429620" & subject_id != "014789935")

################################################################################
### fix date formatting
ff <- ff %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                                        grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                                        T ~ NA_character_))

################################################################################

# create FASTA filename string
ff$FASTAfilename <- paste0(ff$PlateName, ".all.consensus.final.genbank.fasta")

### constants
ff$Subtype <- "B"
ff$Passage <- "Original"


# ### create location from state collection location
# ff <- separate(data = ff, col = SiteName, sep = "\\_", into = c("Site", "StateAbbrev"), fill = "right")
# #ff$State <- state.name[match(ff$StateAbbrev,state.abb)]
# ff <- ff %>% mutate(State = case_when(received_source == "RVTN" ~ Site, 
#                                       received_source == "VIEW" ~ Site, 
#                                       T ~ state.name[match(ff$StateAbbrev,state.abb)]))

# ff <- ff %>% mutate(StateAbbrev = case_when(received_source == "RVTN" ~ state.abb[match(ff$State,state.name)],
#                                             received_source == "VIEW" ~ state.abb[match(ff$State,state.name)], 
#                                       T ~ StateAbbrev))


# create virus name
# hCoV-19/USA/MI-UM-10037140915/2020
# hCoV-19/USA/MA-IVY-ZZX9KKEV/2021
# ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("hRSV/A/USA/", state, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
#                                           received_source == "CDCIVY4" ~ paste0("hRSV/A/USA/", state, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
#                                           received_source == "CDCIVY5" ~ paste0("hRSV/A/USA/", state, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
#                                           #received_source == "RVTN" ~ paste0("hRSV/A/USA/", state, "-RVTN-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
#                                           #received_source == "VIEW" ~ paste0("hRSV/A/USA/", state, "-VIEW-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
#                                           received_source == "IVYIC" ~ paste0("hRSV/A/USA/", state, "-IVYIC-", sample_id, "/", substr(coll_date, 1, 4)),
#                                           received_source == "MDHHS" ~ paste0("hRSV/A/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4)),
#                                           T ~ paste0("hRSV/A/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4))))

ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("hRSV-B-", sample_id),
                                          received_source == "CDCIVY4" ~ paste0("hRSV-B-", sample_id),
                                          received_source == "CDCIVY5" ~ paste0("hRSV-B-", sample_id),
                                          received_source == "CDCIVY6" ~ paste0("hRSV-B-", sample_id),
                                          received_source == "CDCIVY7" ~ paste0("hRSV-B-", sample_id),
                                          #received_source == "RVTN" ~ paste0("hRSV/A/USA/", state, "-RVTN-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
                                          #received_source == "VIEW" ~ paste0("hRSV/A/USA/", state, "-VIEW-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
                                          received_source == "RIGHT" ~ paste0("hRSV-B-", sample_id_lauring),
                                          received_source == "IVYIC" ~ paste0("hRSV-B-", sample_id),
                                          received_source == "MDHHS" ~ paste0("hRSV-B-", sample_id),
                                          T ~ paste0("hRSV-B-", sample_id)))


### constants
ff$Sequence_ID <- ff$VirusName
ff$Collection_date <- paste0(day(ff$coll_date), "-", month(ff$coll_date, label = TRUE), "-", year(ff$coll_date))
ff$Country <- "USA"
ff$Host <- "Human"
ff$Strain <- "B"


################################################################################

### write out VirusName + sample_id crosswalk for use in making 
# .all.consensus.final.gisaid.fasta

ff_crosswalk <- ff %>% select(sample_id, VirusName)

write.csv(ff_crosswalk, paste0(starting_path, "/SEQUENCING/RSV_B/3_ProcessedGenomes/", plate_datef, "_RSVB_", runtech, "_Run_", runnum, "/", plate_datef, "_RSVB_", runtech, "_Run_", runnum, ".forgenbank.meta.csv"), row.names = FALSE, na = "")

## select variables
ff_writeout <- ff %>% select(Sequence_ID, Collection_date, Country, Host, Strain)

ff_writeout <- ff_writeout %>% distinct()

## gisaid upload file name
today <- current_date_string()
gufn <- paste0(today, "_Lauring_genbank_upload_metadata_run_", runnum)

write.table(ff_writeout, paste0(outputLOC, gufn, ".txt"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)

# ## write to excel file (follow format)
# wb <- loadWorkbook(paste0(starting_path, "/SEQUENCING/RSV_A/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/GISAID_UPLOAD_TEMPLATE_2.xlsx"))
# 
# # fill in the submissions tab with built data frame
# writeData(wb, ff_writeout, sheet = "Submissions", startRow = 3, startCol = 1, colNames = FALSE)
# 
# saveWorkbook(wb, paste0(outputLOC, gufn, ".xlsx"), overwrite = TRUE)

