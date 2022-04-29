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
  starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  
} else if (grepl("leighbaker", checking_wd)){
  code_path <- "/Users/leighbaker/Documents/Lauring_Lab/COVID19SequenceCompile/"
  starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"
  
} else {
  
  print("User not recognized.")
  
}


source(paste0(code_path, "pipeline_functions.R"))

################################################################################
### fill in some info

plate_datef <- strsplit(plate_name, "_")[[1]][1] # plate date in YYYYMMDD format
runtech <- strsplit(plate_name, "_")[[1]][3] # nanopore or illumina, will match "PlatePlatform" options
runnum <- strsplit(plate_name, "_")[[1]][5] # number, will match "PlateNumber" options

################################################################################

# create gisaid directory
dir.create(paste0(starting_path, "/SEQUENCING/SARSCOV2/5_GISAID_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum))

################################################################################

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
source(paste0(code_path, "OutsidePipeline/checking_compiled_files.R"))

# set output path for gisaid upload file
# will need to add appropriate folder name at the end of this path
outputLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/5_GISAID_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum, "/")

################################################################################

# read in full compiled pile
finalfileLOC <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary")
final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")

final_file <- filter(final_file, !grepl("Missing Date in Manifest", flag))

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_completeness) >= 90)

# select run of choice
ff <- filter(ff, PlatePlatform == runtech & PlateNumber == runnum)

table(ff$received_source, useNA = "always")
#ff <- filter(ff, received_source != "")
################################################################################
# set up alert to duplicate items

if (any(ff$sample_per_subject > 1)){
  print("STOP: Examine this set of GISAID submissions.")
  stop("There are samples from subject_ids that we've sequenced previously.")
}


#samples_previous <- filter(ff, sample_per_subject > 1) %>% select(subject_id, sample_id, coll_date)
#original_full <- filter(final_file, subject_id %in% unique(samples_previous$subject_id))
### check if the samples are > 90 days apart from one another - then you can let 
### them through.

### uncomment this portion to remove those samples
### to remove these: 
#ff <- filter(ff, sample_per_subject == 1 | !subject_id %in% c("101437962"))
#ff <- filter(ff, sample_per_subject == 1)
#ff <- filter(ff, sample_id != "10041097200")
#ff <- filter(ff, subject_id != "033646963" & subject_id != "025310652")

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
  
} else {
  
  print("User not recognized.")
  
}

# create FASTA filename string
ff$FASTAfilename <- paste0(ff$PlateName, ".all.consensus.final.gisaid.fasta")

### constants
ff$Type <- "betacoronavirus"
ff$Passage <- "Original"

### create location from state collection location
ff <- separate(data = ff, col = SiteName, sep = "\\_", into = c("Site", "StateAbbrev"), fill = "right")
#ff$State <- state.name[match(ff$StateAbbrev,state.abb)]
ff <- ff %>% mutate(State = case_when(received_source == "RVTN" ~ Site, 
                                      T ~ state.name[match(ff$StateAbbrev,state.abb)]))

ff <- ff %>% mutate(StateAbbrev = case_when(received_source == "RVTN" ~ state.abb[match(ff$State,state.name)], 
                                      T ~ StateAbbrev))

ff <- ff %>% mutate(Location = case_when(received_source == "CDCIVY" ~ paste0("North America / USA / ", State), 
                                         received_source == "CDCIVY4" ~ paste0("North America / USA / ", State), 
                                         received_source == "RVTN" ~ paste0("North America / USA / ", State), 
                                         T ~ "North America / USA / Michigan"))

# create virus name
# hCoV-19/USA/MI-UM-10037140915/2020
# hCoV-19/USA/MA-IVY-ZZX9KKEV/2021
ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("hCoV-19/USA/", StateAbbrev, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
                                          received_source == "CDCIVY4" ~ paste0("hCoV-19/USA/", StateAbbrev, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)),
                                          received_source == "RVTN" ~ paste0("hCoV-19/USA/", StateAbbrev, "-RVTN-", sample_id, "/", substr(coll_date, 1, 4)),
                                          T ~ paste0("hCoV-19/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4))))


### constants
ff$AdditionalLoc <- ""
ff$Host <- "Human"
ff$AdditionalHost <- ""

ff <- ff %>% mutate(SamplingStrategy = case_when(sample_per_subject > 1 ~ "Warning", 
                                                 received_source %in% c("CDCIVY", "MHOME") ~ "", 
                                                 grepl("PUI", flag) ~ "", 
                                                 received_source == "RVTN" ~ "Research",
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
ff$SequencingTechnology <- ifelse(ff$PlatePlatform == "Nanopore", "Oxford Nanopore", 
                                  ifelse(ff$PlatePlatform == "Illumina", "Illumina MiSeq", "Unknown"))

unknown_tech <- filter(ff, SequencingTechnology == "Unknown")

if (nrow(unknown_tech) != 0){
  stop("Check Sequencing Technology options.")
}

### Assembly Method
ff$AssemblyMethod <- ifelse(ff$PlatePlatform == "Nanopore", "ARTIC Network pipeline", 
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
                                          received_source == "RVTN" ~ "Vanderbilt University Medical Center",
                                          T ~ "University of Michigan Clinical Microbiology Laboratory"), 
                    originlabaddress = case_when(received_source == "CDCIVY" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "CDCIVY4" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA",
                                                 received_source == "RVTN" ~ "Medical Center North CC303, 1161 21st Ave. S., Nashville, TN, USA",
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
  
} else {
  
  print("User not recognized.")
  
}

ff$comment <- ""
ff$commenticon <- ""

################################################################################

### write out VirusName + sample_id crosswalk for use in making 
# .all.consensus.final.gisaid.fasta

ff_crosswalk <- ff %>% select(sample_id, VirusName)
write.csv(ff_crosswalk, paste0(starting_path, "/SEQUENCING/SARSCOV2/3_ProcessedGenomes/", plate_datef, "_SC2_", runtech, "_Run_", runnum, "/", plate_datef, "_SC2_", runtech, "_Run_", runnum, ".forgisaid.meta.csv"), row.names = FALSE, na = "")

## select variables
ff_writeout <- ff %>% select(Submitter, FASTAfilename, VirusName,Type, Passage,  coll_date, Location, 
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
wb <- loadWorkbook(paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/GISAID_UPLOAD_TEMPLATE2.xlsx"))

# fill in the submissions tab with built data frame
writeData(wb, ff_writeout, sheet = "Submissions", startRow = 3, startCol = 1, colNames = FALSE)

saveWorkbook(wb, paste0(outputLOC, gufn, ".xlsx"), overwrite = TRUE)

