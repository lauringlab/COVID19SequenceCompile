################################################################################
#               Complete Run - COVID-19 Genetic Data Compilation               #
#                         Last Updated: 07/01/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
batch_path <- "C:/Users/juliegil/Documents/UofM_Work/Lab_Organization/AlertCode"

################################################################################

source(paste0(code_path, "pipeline_functions.R"))

################################################################################

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))

source(paste0(code_path, "pangolin_code.R"))

source(paste0(code_path, "nextclade_code.R"))

source(paste0(code_path, "gisaid_code.R"))

source(paste0(code_path, "compile_components_code.R"))

shell.exec(paste0(batch_path, "/sample_count_run.bat"))

