################################################################################
#               GISAID File Upload Format Creation for Influenza               #
#                            Created: 11/19/2021                               #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(openxlsx)

################################################################################
# just need some of these functions

code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
source(paste0(code_path, "pipeline_functions.R"))

# set starting path
starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

################################################################################
### fill in some info manually

plate_datef <- "20211111" # plate date in YYYYMMDD format
runtech <- "Nanopore" # nanopore or illumina, will match "PlatePlatform" options
runnum <- "3" # number, will match "PlateNumber" options

################################################################################

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
code_path2 <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/InfluenzaACode/"
source(paste0(code_path2, "OutsidePipeline/comparing_full_secret_influenza.R"))

# set output path for gisaid upload file
# will need to add appropriate folder name at the end of this path
outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/5_GISAID_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum, "/")

################################################################################

# read in full compiled pile
finalfileLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary")
final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_HA_completeness) >= 90)

# select run of choice
ff <- filter(ff, PlatePlatform == runtech & PlateNumber == runnum)

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
#ff <- filter(ff, sample_per_subject == 1)
#ff <- filter(ff, sample_per_subject == 1 | sample_id == "10042234896")
#ff <- filter(ff, sample_id != "10041097200")

################################################################################
### fix date formatting
ff <- ff %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

################################################################################

### strain naming: A/Michigan/UOM[sample_id]/2021
### {flu type A/B} / {collection state} / {UOM}{sample_id} / {collection year}
ff$StrainName <- ifelse(ff$received_source %in% c("CBR", "UHS"), paste0("A/Michigan/UOM", ff$sample_id, "/", year(ff$coll_date)), "CHECK")

if (any(ff$StrainName == "CHECK")){
  stop("Unexpected received source!")
}

ff$Passage <- "ORI"
ff$Type <- "A"
ff$SubtypeH <- gsub("H", "", ff$nextclade_HA_type)
ff$SubtypeN <- ifelse(ff$nextclade_HA_type == "H3", "2", #H3N2
                      ifelse(ff$nextclade_HA_type == "H1", "1", "CHECK")) #H1N1

if(any(ff$SubtypeN == "CHECK")){
  stop("Issue with N subtype assignment")
}

ff$CollectionDate <- ff$coll_date
ff$Location <- "North America"
ff$Host <- "Human"

#University of Michigan Clinical Microbiology Laboratory
#2800 Plymouth Rd, Ann Arbor, MI, USA

ff_gisaid <- ff %>% select(StrainName, Passage, Type, SubtypeH, SubtypeN, 
                           CollectionDate, Location, Host)

## single upload: A/Michigan/UOM10042526240/2021 (2021-11-21)