################################################################################
#                    GISAID File Upload Format Creation                        #
#                         Last Updated: 05/25/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(openxlsx)

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

# read in full compiled pile
finalfileLOC <- paste0(starting_path, "SequenceSampleMetadata/FinalSummary")
final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_completeness) >= 90)

# select run of choice
ff <- filter(ff, PlatePlatform == "" & PlateNumber == "")

# enter GISAID username here
ff$Submitter <- ""

# create FASTA filename string
ff$FASTAfilename <- paste0(gsub("-", "", ff$PlateDate), "_", ff$PlateName, "_", ff$PlateNumber, ".all.consensus.final.gisaid.fasta")

# create virus name
# hCoV-19/USA/MI-UM-10037140915/2020

ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("hCoV-19/USA/CDC-IVY-", sample_id, "/", substr(ff$coll_date, 1, 4)), 
                                          T ~ paste0("hCoV-19/USA/MI-UM-", sample_id, "/", substr(ff$coll_date, 1, 4))))

### constants
ff$Type <- "betacoronavirus"
ff$Passage <- "Original"

### create location from state collection location
ff <- separate(data = ff, col = SiteName, sep = "\\_", into = c("Site", "StateAbbrev"))
ff$State <- state.name[match(ff$StateAbbrev,state.abb)]

ff <- ff %>% mutate(Location = case_when(received_source == "CDCIVY" ~ paste0("North America / USA / ", State), 
                                         T ~ "North America / USA / Michigan"))

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
