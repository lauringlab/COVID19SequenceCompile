library(tidyverse)
library(openxlsx)

# get plate map from 
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\2_PlateMaps
# and drop it in
# Dropbox (University of Michigan)\MED-LauringLab\SEQUENCING\SARSCOV2\4_SequenceSampleMetadata\PlateMaps

paste0(starting_path, "SEQUENCING/SARSCOV2/2_PlateMaps/", plate_name, ".xlsx")

file.copy(from = paste0(starting_path, "SEQUENCING/SARSCOV2/2_PlateMaps/", plate_name, ".xlsx"),
            to = paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/PlateMaps/", plate_name, ".xlsx"))

