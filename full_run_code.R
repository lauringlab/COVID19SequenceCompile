################################################################################
#               Complete Run - COVID-19 Genetic Data Compilation               #
#                         Last Updated: 07/01/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
batch_path <- "C:/Users/juliegil/Documents/UofM_Work/Lauring_Lab/Lab_Organization/AlertCode"
influenza_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/InfluenzaACode/"

################################################################################

options(scipen=999)
source(paste0(code_path, "pipeline_functions.R"))

################################################################################

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "pangolin_code.R"))

source(paste0(code_path, "nextclade_code.R"))

source(paste0(code_path, "gisaid_code.R"))

source(paste0(code_path, "compile_components_code.R"))

## only run for updating Slack channel (underworld) with total sample count
shell.exec(paste0(batch_path, "/sample_count_run.bat"))

################################################################################

source(pate0(influenza_path, "influenza_manifestcode.R"))

source(pate0(influenza_path, "influenza_platemapcode.R"))

################################################################################


