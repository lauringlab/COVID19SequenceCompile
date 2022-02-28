################################################################################
#               Complete Run - COVID-19 Genetic Data Compilation               #
#                         Last Updated: 02/03/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

#starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
#code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
#batch_path <- "C:/Users/juliegil/Documents/UofM_Work/Lauring_Lab/Lab_Organization/AlertCode"
batch_path <- "/Users/juliegil/Documents/LauringLab_Code/AlertCode"
#influenza_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/InfluenzaACode/"
influenza_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/InfluenzaACode/"

################################################################################

options(scipen=999)
source(paste0(code_path, "pipeline_functions.R"))

################################################################################

# sars-cov-2 plate
plate_name <- "20220210_SC2_Nanopore_Run_120"

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

source(paste0(code_path, "compile_components_code.R"))

source(paste0(code_path, "OutsidePipeline/subset_compiled_for_fasta.R"))

################################################################################

# checking_sampleids.py step goes here
# may need keyboard shortcut = Ctrl+Shift+F10 to clear environment/loaded packages
#system(paste0("c:/users/juliegil/appdata/local/programs/python/python38/python.exe ", code_path, "OutsidePipeline/ProcessingFASTA/checking_sampleids.py ", plate_name))
#system("pip3 install pandas")
#system(paste0("python3 ", code_path, "OutsidePipeline/ProcessingFASTA/checking_sampleids.py ", plate_name))


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

# move into the platemap folder of the run you're on.

# Activate the pangolin environment 
# The command is:
# conda activate pangolin

# Run pangolin
# The command is: 
# pangolin {plateMapName}.all.consensus.renamed.full.fasta

# This will output a file called lineage_report.csv 

# Rename this file to plateMapName_pangolin.csv
# Copy this renamed file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/pangolin]
source(paste0(code_path, "OutsidePipeline/moving_pangolin_output.R"))

################################################################################
#                                NEXTCLADE                                     #
################################################################################

# Steps:

# Navigate to https://clades.nextstrain.org/

# In the SARS-CoV-2 box in the lower right-hand corner, ensure the "From File" 
# tab is selected and either drag & drop or select the 
# plateMapName.all.consensus.renamed.full.fasta from the 
# [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] folder you've been working from

# Once complete, in the top bar of the processing screen, click the arrow next 
# to "Export to CSV" and select "Export to TSV"

# The file will go to the Downloads file of your computer as nextclade.tsv - 
# Copy the file to [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] and 
# rename it to plateMapName_nextclade.tsv

# Copy that re-named file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/nextclade]
source(paste0(code_path, "OutsidePipeline/moving_nextclade_output.R"))


### nextclade CLI download information: https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html

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
#shell.exec(paste0("sh ", batch_path, "/sample_count_run_mac.sh"))
#system2(paste0("sh ", batch_path, "/sample_count_run_mac.sh"))

### next, do gisaid upload steps

################################################################################
#                                 ROUND 3                                      #
################################################################################
# Run this after downloading new gisaid metadata from gisaid.org

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


