################################################################################
#               Complete Run - COVID-19 Genetic Data Compilation               #
#                         Last Updated: 11/19/2021                             #
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

# sars-cov-2 plate
plate_name <- "20220127_SC2_Nanopore_Run_108"

################################################################################
#                                 ROUND 1                                      #
################################################################################

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "compile_components_code.R"))

source(paste0(code_path, "OutsidePipeline/subset_compiled_for_fasta.R"))

# checking_sampleids.py step goes here

################################################################################
#                                 ROUND 2                                      #
################################################################################

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "pangolin_code.R"))

source(paste0(code_path, "nextclade_code.R"))

source(paste0(code_path, "compile_components_code.R"))

## only run for updating Slack channel (underworld) with total sample count
shell.exec(paste0(batch_path, "/sample_count_run.bat"))

################################################################################
#                                 ROUND 3                                      #
################################################################################

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "pangolin_code.R"))

source(paste0(code_path, "nextclade_code.R"))

source(paste0(code_path, "gisaid_code.R"))

source(paste0(code_path, "compile_components_code.R"))

################################################################################

################################################################################
# INFLUENZA COMPONENT
################################################################################

source(paste0(influenza_path, "influenza_manifestcode.R"))

source(paste0(influenza_path, "influenza_platemapcode.R"))

source(paste0(influenza_path, "influenza_nextcladecode.R"))

### will need a piece to incorporate gisaid returns

### need to compile everything
source(paste0(influenza_path, "influenza_compilecomponentscode.R"))

################################################################################


