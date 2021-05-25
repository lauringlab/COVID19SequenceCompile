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


