### Code to generate RedCap upload file for RSVA - IVY

library(tidyverse)
library(lubridate)


# read in full compiled rsv_a file

rsv_a <- read.csv("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv", colClasses = "character")
