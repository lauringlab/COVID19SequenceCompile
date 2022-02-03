library(tidyverse)
library(openxlsx)

# get plate map from 
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\2_PlateMaps
# and drop it in
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\4_SequenceSampleMetadata\PlateMaps

file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/2_PlateMaps/", plate_name, ".xlsx"),
            to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name, ".xlsx"))


# then, read in that excel file that we moved to the inner pipeline set
file_in <- read.xlsx(paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name, ".xlsx"), sheet = 1)

# remove the well position column
drop_columns <- c("Well.Position")
file_in <- file_in[!names(file_in) %in% drop_columns]

# Rename first 7 columns
colnames(file_in)[1] <- "Processing Plate"
colnames(file_in)[2] <- "Slot"
colnames(file_in)[3] <- "Sample ACCN"
colnames(file_in)[4] <- "Sample MRN"
colnames(file_in)[5] <- "Sample Order#"
colnames(file_in)[6] <- "Barcode"
colnames(file_in)[7] <- "Source"

# write the excel file back out
write.xlsx(file_in, paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name, ".xlsx"), overwrite = TRUE)
