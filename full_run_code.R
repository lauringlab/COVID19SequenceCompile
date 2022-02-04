################################################################################
#               Complete Run - COVID-19 Genetic Data Compilation               #
#                         Last Updated: 02/03/2021                             #
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
plate_name <- "20220128_SC2_Nanopore_Run_112"

################################################################################
#                                 ROUND 1                                      #
################################################################################
# Explanation: Moves proper plate map from load location to use location, compiles
# manifests and plate maps, merges them together, and checks that the plate we're 
# running has a manifest entry for each sample
################################################################################

source(paste0(code_path, "OutsidePipeline/moving_plate_map_files.R"))

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open
# you'll get this error if this is the case: 
# Error in file(con, "r") : invalid 'description' argument
# In addition: Warning message:
# In unzip(xlsxFile, exdir = xmlDir) : error 1 in extracting from zip file

source(paste0(code_path, "compile_components_code.R"))

source(paste0(code_path, "OutsidePipeline/subset_compiled_for_fasta.R"))
## Common error here: 
#Error in file(file, ifelse(append, "a", "w")) : 
#  cannot open the connection
#In addition: Warning message:
#  In file(file, ifelse(append, "a", "w")) :
#  cannot open file 'C:/Users/juliegil/Dropbox (University of Michigan)/
# MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes/20220128_SC2_Nanopore_Run_112/
# 20220128_SC2_Nanopore_Run_112.meta.csv': No such file or directory

# usually means that the plate map name, and the folder created in 3_ProcessedGenomes
# for a particular run, don't match

################################################################################

# checking_sampleids.py step goes here
# may need keyboard shortcut = Ctrl+Shift+F10 to clear environment/loaded packages
system(paste0("c:/users/juliegil/appdata/local/programs/python/python38/python.exe ", code_path, "OutsidePipeline/ProcessingFASTA/checking_sampleids.py ", plate_name))

################################################################################
#                                PANGOLIN                                      #
################################################################################

# Instructions for installing the command line version of Pangolin can be found 
# here: https://cov-lineages.org/pangolin_docs/installation.html

# Note: This can only be completed using Linux or OS operating systems. (The 
# environment.yml file to install pangolin requires minimap2 and gofasta libraries, 
# both of which are only built for those systems, not for Windows) If you only 
# have a Windows operating system, you can install a Windows Subsystem for Linux 
# (directions here: https://docs.microsoft.com/en-us/windows/wsl/install-win10) 
# which will allow you to proceed with the pangolin installation process.

# Steps:
  
# In command line, navigate to [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName]
# For me, this command is: 
# cd /mnt/c/Users/juliegil/'Dropbox (University of Michigan)'/MED-LauringLab/SEQUENCING/SARSCOV2/3_ProcessedGenomes

# Activate the pangolin environment 
# The command is:
# conda activate pangolin

# Run pangolin
# The command is: 
# pangolin {plateMapName}.all.consensus.renamed.full.fasta

# This will output a file called lineage_report.csv 

# Rename this file to plateMapName_pangolin.csv

# Copy this renamed file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/pangolin]

################################################################################
#                                NEXTCLADE                                     #
################################################################################




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


