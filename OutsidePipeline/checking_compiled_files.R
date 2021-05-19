################################################################################
#          Check of Cumulative Dataset for COVID-19 Genetic Sampling           #
#                         Last Updated: 05/18/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
compiled_loc <- paste0(starting_path, "SequenceSampleMetadata/FinalSummary")

# essentially, we just want a small code snippet to compare the "secret" version 
# of the complete file, to the "accessible" version of the complete file

accessible_version <- read.csv(paste0(compiled_loc, "/full_compiled_data.csv"), colClasses = "character")
secret_version <- read.csv(paste0(compiled_loc, "/secret/full_compiled_data.csv"), colClasses = "character")

## this is the check, comparing two data frames
x <- all_equal(accessible_version, secret_version)

## here is the information: 
if (x){
  print("The Secret Version and the Accessible Version of full_compiled_data.csv are the same.")
} else {
  print("The Secret Version and the Accessible Version of full_compiled_data.csv are not the same.")
  print("Differences:")
  print(x)
}


