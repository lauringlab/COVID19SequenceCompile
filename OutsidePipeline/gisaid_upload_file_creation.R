################################################################################
#                    GISAID File Upload Format Creation                        #
#                         Last Updated: 05/27/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(openxlsx)

################################################################################
# just need some of these functions

code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
source(paste0(code_path, "pipeline_functions.R"))

################################################################################
### fill in some info manually

plate_datef <- "20210524" # plate date in YYYYMMDD format
runtech <- "Nanopore" # nanopore or illumina, will match "PlatePlatform" options
runnum <- "27" # number, will match "PlateNumber" options

################################################################################

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
source("C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/OutsidePipeline/checking_compiled_files.R")

# set starting path
starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

# set output path for gisaid upload file
# will need to add appropriate folder name at the end of this path
outputLOC <- paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/GISAID_Uploads/upload_", plate_datef, "_", tolower(runtech), "_run_", runnum, "/")

################################################################################

# read in full compiled pile
finalfileLOC <- paste0(starting_path, "SequenceSampleMetadata/FinalSummary")
final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_completeness) >= 90)

# select run of choice
ff <- filter(ff, PlatePlatform == runtech & PlateNumber == runnum)

#table(final_file$PlatePlatform, final_file$PlateNumber, useNA = "always")
################################################################################
# need to compare to complete GISAID file, to avoid submitting duplicate sequences

################################################################################

# enter GISAID username here
ff$Submitter <- "juliegil"

# create FASTA filename string
ff$FASTAfilename <- paste0(ff$PlateName, ".all.consensus.final.gisaid.fasta")

### constants
ff$Type <- "betacoronavirus"
ff$Passage <- "Original"

### create location from state collection location
ff <- separate(data = ff, col = SiteName, sep = "\\_", into = c("Site", "StateAbbrev"))
ff$State <- state.name[match(ff$StateAbbrev,state.abb)]

ff <- ff %>% mutate(Location = case_when(received_source == "CDCIVY" ~ paste0("North America / USA / ", State), 
                                         T ~ "North America / USA / Michigan"))

# create virus name
# hCoV-19/USA/MI-UM-10037140915/2020
# hCoV-19/USA/MA-IVY-ZZX9KKEV/2021
ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("hCoV-19/USA/", StateAbbrev, "-IVY-", sample_id, "/", substr(coll_date, 1, 4)), 
                                          T ~ paste0("hCoV-19/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4))))


### constants
ff$AdditionalLoc <- ""
ff$Host <- "Human"
ff$AdditionalHost <- ""
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
                                          T ~ "University of Michigan Clinical Microbiology Laboratory"), 
                    originlabaddress = case_when(received_source == "CDCIVY" ~ "Medical Center North D7240, 1161 21st Ave. S., Nashville, TN, USA", 
                                                  T ~ "2800 Plymouth Rd, Ann Arbor, MI, USA"))

ff$originlabsampleid <- ""

### submitting Lab
ff$submitlab <- "Lauring Lab, University of Michigan, Department of Microbiology and Immunology"
ff$submitlabaddress <- "1150 W. Medical Center Dr. MSRB1, Ann Arbor, MI, USA"
ff$submitlabsampleid <- ""

### Authors
ff$authors <- "Gilbert"

ff$comment <- ""
ff$commenticon <- ""

################################################################################

### write out VirusName + sample_id crosswalk for use in making 
# .all.consensus.final.gisaid.fasta

ff_crosswalk <- ff %>% select(sample_id, VirusName)
write.csv(ff_crosswalk, paste0(starting_path, "/ProcessedGenomes/", plate_datef, "_", runtech, "_Run_", runnum, "/", plate_datef, "_", runtech, "_Run_", runnum, ".forgisaid.meta.csv"), row.names = FALSE, na = "")

## select variables
ff_writeout <- ff %>% select(Submitter, FASTAfilename, VirusName,Type, Passage,  coll_date, Location, 
                             AdditionalLoc, Host, AdditionalHost, Gender, Age, Status, 
                             SpecimenSource, Outbreak, lastVaccinated, Treatment, SequencingTechnology, 
                             AssemblyMethod, Coverage, originlab, originlabaddress, originlabsampleid, 
                             submitlab, submitlabaddress, submitlabsampleid, authors, 
                             comment, commenticon)

ff_writeout <- ff_writeout %>% distinct()

## gisaid upload file name
today <- current_date_string()
gufn <- paste0(today, "_Lauring_gisaid_upload_metadata_run_", runnum) 

## write to excel file (follow format)
wb <- loadWorkbook(paste0(starting_path, "/SequenceSampleMetadata/SequenceOutcomes/gisaid/GISAID_UPLOAD_TEMPLATE.xlsx"))

# fill in the submissions tab with built data frame
writeData(wb, ff_writeout, sheet = "Submissions", startRow = 3, startCol = 1, colNames = FALSE)

saveWorkbook(wb, paste0(outputLOC, gufn, ".xlsx"), overwrite = TRUE)
