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
  
  #code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
  code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
  starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  #code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
  #starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  
} else if (grepl("leighbaker", checking_wd)){
  code_path <- "/Users/leighbaker/Documents/Lauring_Lab/COVID19SequenceCompile/"
  starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"

} else if (grepl("leighbak", checking_wd)){
  
  starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"
  code_path <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/"
  batch_path <- "/Users/leighbak/Documents/Lauring_Lab/AlertCode"
    
} else {
  
  print("User not recognized.")
  
}


source(paste0(code_path, "pipeline_functions.R"))

################################################################################
### fill in some info

#fill in the plate name below if running this code seperate and not after "full_run_code.R"
#plate_name <- "20230323_RSVB_Illumina_Run_1"

plate_datef <- strsplit(plate_name, "_")[[1]][1] # plate date in YYYYMMDD format
runtech <- strsplit(plate_name, "_")[[1]][3] # nanopore or illumina, will match "PlatePlatform" options
runnum <- strsplit(plate_name, "_")[[1]][5] # number, will match "PlateNumber" options

################################################################################

# create gisaid directory
dir.create(paste0(starting_path, "SEQUENCING/RSV_B/5_GISAID_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum))

################################################################################

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
source(paste0(code_path, "OutsidePipeline/checking_compiled_files.R"))

# set output path for gisaid upload file
# will need to add appropriate folder name at the end of this path
outputLOC <- paste0(starting_path, "SEQUENCING/RSV_B/5_GISAID_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum, "/")

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
  print("STOP: Examine this set of GISAID submissions.")
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

# enter GISAID username here 
ff$Submitter <- 
  if (grepl("juliegil", checking_wd)){
  
  Submitter <- "juliegil"
  
} else if (grepl("leighbaker", checking_wd)){
  Submitter <- "Leighbaker"
  
} else if (grepl("leighbak", checking_wd)){
  Submitter <- "Leighbaker"
  
} else {
  
  print("User not recognized.")
  
}

# create FASTA filename string
ff$FASTAfilename <- paste0(ff$PlateName, ".all.consensus.final.gisaid.fasta")

### constants
ff$Subtype <- "B"
ff$Passage <- "Original"


if (any(grepl("IVY", ff$received_source))){
  
  ff$source_state_code <- substr(ff$subject_id, 3, 4)
  cdcivy_manifest_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/CDCIVY")
  cdc_sites <- read.csv(paste0(cdcivy_manifest_fp, "/Keys/CDC_SiteCodebook.csv"), colClasses = "character") %>% separate(SiteCode, into = c("site", "state"), sep = "_")
  site_bit <- cdc_sites %>% select(Number, state)
  # add leading zeros to site
  site_bit <- site_bit %>% mutate(Number = case_when(nchar(Number) == 1 ~ paste0("0", Number), 
                                                     T ~ Number))
  #colnames(site_bit) <- c("Number", "state")
  ff <- merge(ff, site_bit, by.x = c("source_state_code"), by.y = c("Number"), all.x = TRUE)
  ff$state <- ifelse(!grepl("IVY", ff$received_source), "", ff$state)
  
} else {
  state <- ""
}


# ### create location from state collection location
# ff <- separate(data = ff, col = SiteName, sep = "\\_", into = c("Site", "StateAbbrev"), fill = "right")
# #ff$State <- state.name[match(ff$StateAbbrev,state.abb)]
# ff <- ff %>% mutate(State = case_when(received_source == "RVTN" ~ Site, 
#                                       received_source == "VIEW" ~ Site, 
#                                       T ~ state.name[match(ff$StateAbbrev,state.abb)]))

# ff <- ff %>% mutate(StateAbbrev = case_when(received_source == "RVTN" ~ state.abb[match(ff$State,state.name)],
#                                             received_source == "VIEW" ~ state.abb[match(ff$State,state.name)], 
#                                       T ~ StateAbbrev))

ff <- ff %>% mutate(Location = case_when(received_source == "CDCIVY" ~ paste0("North America / USA / ", state.name[match(ff$state,state.abb)]), 
                                         received_source == "CDCIVY4" ~ paste0("North America / USA / ", state.name[match(ff$state,state.abb)]), 
                                         received_source == "CDCIVY5" ~ paste0("North America / USA / ", state.name[match(ff$state,state.abb)]), 
                                         received_source == "IVYIC" ~ paste0("North America / USA / ", state.name[match(ff$state,state.abb)]),
                                         received_source == "RVTN" ~ paste0("North America / USA / ", state.name[match(ff$state,state.abb)]), 
                                         received_source == "VIEW" ~ paste0("North America / USA / ", state.name[match(ff$state,state.abb)]), 
                                         T ~ "North America / USA / Michigan"))

# create virus name
# hCoV-19/USA/MI-UM-10037140915/2020
# hCoV-19/USA/MA-IVY-ZZX9KKEV/2021
ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("hRSV/B/USA/", state, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
                                          received_source == "CDCIVY4" ~ paste0("hRSV/B/USA/", state, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
                                          received_source == "CDCIVY5" ~ paste0("hRSV/B/USA/", state, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
                                          #received_source == "RVTN" ~ paste0("hRSV/A/USA/", state, "-RVTN-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
                                          #received_source == "VIEW" ~ paste0("hRSV/A/USA/", state, "-VIEW-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
                                          received_source == "IVYIC" ~ paste0("hRSV/B/USA/", state, "-IVYIC-", sample_id, "/", substr(coll_date, 1, 4)),
                                          received_source == "MDHHS" ~ paste0("hRSV/B/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4)),
                                          T ~ paste0("hRSV/B/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4))))


### constants
ff$AdditionalLoc <- ""
ff$Host <- "Human"
ff$AdditionalHost <- ""

ff <- ff %>% mutate(SamplingStrategy = case_when(sample_per_subject > 1 ~ "Warning", 
                                                 received_source %in% c("CDCIVY", "CDCIVY4", "CDCIVY5", "MHOME") ~ "", 
                                                 grepl("PUI", flag) ~ "", 
                                                 received_source == "RVTN" ~ "Research",
                                                 received_source == "VIEW" ~ "Research",
                                                 received_source == "IVYIC" ~ "Serial sampling",
                                                 T ~ "Baseline surveillance"))

if(any(ff$SamplingStrategy == "Warning")){
  print("Look at this, apply logic if necessary")
  ## For someone positive > 90 days apart it's tricky. Strictly speaking, 
  ## they would be considered possible reinfection and not longitudinal (and 
  ## therefore surveillance). However, they could be prolonged shedder. One way 
  ## around this would be to look at the lineage in the duplicates, if different, 
  ## then it is reinfection and surveillance for both. Put another way, the 
  ## filter for duplicates would be check dates (>90 days) and check lineage 
  ## (different) in order for it to be considered surveillance.
  
  #ff$SamplingStrategy <- ifelse(ff$SamplingStrategy == "Warning", "Baseline surveillance", ff$SamplingStrategy)
}

#table(ff$SamplingStrategy)

ff$Gender <- "unknown"
ff$Age <- "unknown"
ff$Status <- "unknown"
ff$SpecimenSource <- "unknown"
ff$Outbreak <- ""
ff$lastVaccinated <- ""
ff$Treatment <- ""


# Oxford Nanopore, Illumina MiSeq
ff$SequencingTechnology <- ifelse(ff$PlatePlatform == "Nanopore", "Oxford Nanopore Midnight", 
                                  ifelse(ff$PlatePlatform == "Illumina", "Illumina NextSeq 1000", "Unknown"))

unknown_tech <- filter(ff, SequencingTechnology == "Unknown")

if (nrow(unknown_tech) != 0){
  stop("Check Sequencing Technology options.")
}

### Assembly Method
ff$AssemblyMethod <- ifelse(ff$PlatePlatform == "Nanopore", "ARTIC Network pipeline Midnight", 
                            ifelse(ff$PlatePlatform == "Illumina", "BWA-MEM, iVar", "Unknown"))

unknown_assembly <- filter(ff, AssemblyMethod == "Unknown")

if (nrow(unknown_assembly) != 0){
  stop("Check Assembly Method options.")
}

### Coverage
ff$Coverage <- ""

### Originating Lab
ff <- ff %>% mutate(originlab = case_when(received_source == "CDCIVY" ~ "IVY3 Central Lab, Vanderbilt University Medical Center", 
                                          received_source == "CDCIVY4" ~ "IVY4 Central Lab, Vanderbilt University Medical Center",
                                          received_source == "CDCIVY5" ~ "IVY5 Central Lab, Vanderbilt University Medical Center",
                                          received_source == "RVTN" ~ "Vanderbilt University Medical Center",
                                          received_source == "VIEW" ~ "Vanderbilt University Medical Center",
                                          received_source == "IVYIC" ~ "IVY4 Central Lab, Vanderbilt University Medical Center",
                                          received_source == "MDHHS" ~ "Michigan Department of Health and Human Services, Bureau of Laboratories",
                                          received_source == "TRINITY" ~ "Warde Medical Laboratory",
                                          received_source == "ASC" ~ "TMCIDR Lab",
                                          received_source == "ASJ" ~ "TMCIDR Lab",
                                          received_source == "HFHS" ~ "Henry Ford Health Microbiology Laboratory",
                                          T ~ "University of Michigan Clinical Microbiology Laboratory"), 
                    originlabaddress = case_when(received_source == "CDCIVY" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "CDCIVY4" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "CDCIVY5" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "RVTN" ~ "Medical Center North CC303, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "VIEW" ~ "Medical Center North CC303, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "IVYIC" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "MDHHS" ~ "3350 N Martin Luther King Jr Blvd",
                                                 received_source == "TRINITY" ~ "300 West Textile Rd, Ann Arbor, MI 48108",
                                                 received_source == "ASC" ~ "19251 Mack Avenue Suite 575, Grosse Pointe Woods, MI 48236",
                                                 received_source == "ASJ" ~ "19251 Mack Avenue Suite 575, Grosse Pointe Woods, MI 48236",
                                                 received_source == "HFHS" ~ "2799 West Grand Blvd., 6065 E&R Building, Detroit, MI 48202",
                                                  T ~ "2800 Plymouth Rd, Ann Arbor, MI, USA"))

ff$originlabsampleid <- ""

### submitting Lab
ff$submitlab <- "Lauring Lab, University of Michigan, Department of Microbiology and Immunology"
ff$submitlabaddress <- "1137 Catherine Street, Ann Arbor, MI, USA"
ff$submitlabsampleid <- ""

### Authors
ff$authors <-
if (grepl("juliegil", checking_wd)){
  
  authors <- "Gilbert"
  
} else if (grepl("leighbaker", checking_wd)){
  authors <- "Baker"
  
} else if (grepl("leighbak", checking_wd)){
  authors <- "Baker" 
  
} else {
  
  print("User not recognized.")
  
}

ff$comment <- ""
ff$commenticon <- ""

################################################################################

### write out VirusName + sample_id crosswalk for use in making 
# .all.consensus.final.gisaid.fasta

ff_crosswalk <- ff %>% select(sample_id, VirusName)

write.csv(ff_crosswalk, paste0(starting_path, "/SEQUENCING/RSV_B/3_ProcessedGenomes/", plate_datef, "_RSVB_", runtech, "_Run_", runnum, "/", plate_datef, "_RSVB_", runtech, "_Run_", runnum, ".forgisaid.meta.csv"), row.names = FALSE, na = "")

## select variables
ff_writeout <- ff %>% select(Submitter, FASTAfilename, VirusName,Subtype, Passage,  coll_date, Location, 
                             AdditionalLoc, Host, AdditionalHost, SamplingStrategy, Gender, Age, Status, 
                             SpecimenSource, Outbreak, lastVaccinated, Treatment, SequencingTechnology, 
                             AssemblyMethod, Coverage, originlab, originlabaddress, originlabsampleid, 
                             submitlab, submitlabaddress, submitlabsampleid, authors, 
                             comment, commenticon)

ff_writeout <- ff_writeout %>% distinct()

## gisaid upload file name
today <- current_date_string()
gufn <- paste0(today, "_Lauring_gisaid_upload_metadata_run_", runnum) 

## write to excel file (follow format)
wb <- loadWorkbook(paste0(starting_path, "/SEQUENCING/RSV_A/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/GISAID_UPLOAD_TEMPLATE_2.xlsx"))

# fill in the submissions tab with built data frame
writeData(wb, ff_writeout, sheet = "Submissions", startRow = 3, startCol = 1, colNames = FALSE)

saveWorkbook(wb, paste0(outputLOC, gufn, ".xlsx"), overwrite = TRUE)

