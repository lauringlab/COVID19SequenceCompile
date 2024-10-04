################################################################################
#          Check of Cumulative Dataset for COVID-19 Genetic Sampling           #
#                         Last Updated: 05/18/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(lubridate)

################################################################################

checking_wd <- getwd()
if (grepl("juliegil", checking_wd)){
  
  starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

} else if (grepl("leighbaker", checking_wd)){
  
  starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"
  
} else if (grepl("leighbak", checking_wd)){
  
  starting_path <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/"  
 
} else if (grepl("chbla", checking_wd)){
  
  starting_path <- "/Users/chbla/University of Michigan Dropbox/MED-LauringLab/" 

} else if (grepl("chbl", checking_wd)){
  
  starting_path <- "/Users/chbl/University of Michigan Dropbox/MED-LauringLab/"  
  
} else {
  
  print("User not recognized.")
  
}


compiled_loc <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/FinalSummary")

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


