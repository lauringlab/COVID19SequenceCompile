################################################################################
#               Complete Run - COVID-19 Genetic Data Compilation               #
#                         Last Updated: 07/20/2023                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

checking_wd <- getwd()

if (grepl("juliegil", checking_wd)){
  starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  #starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
  code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
  #code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
  batch_path <- "C:/Users/juliegil/Documents/UofM_Work/Lauring_Lab/Lab_Organization/AlertCode"
  #batch_path <- "/Users/juliegil/Documents/LauringLab_Code/AlertCode"
  influenza_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/InfluenzaACode/"
  #influenza_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/InfluenzaACode/"
  
} else if (grepl("leighbaker", checking_wd)){
  
  starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"
  code_path <- "/Users/leighbaker/Documents/Lauring_Lab/COVID19SequenceCompile/"
  batch_path <- "/Users/leighbaker/Documents/Lauring_Lab/AlertCode"
  
} else if (grepl("leighbak", checking_wd)){
  
  starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"
  code_path <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/"
  batch_path <- "/Users/leighbak/Documents/Lauring_Lab/AlertCode"
  influenza_path <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/InfluenzaACode/"

} else if (grepl("chbla", checking_wd)){
  
  starting_path <- "C:/Users/chbla/Dropbox (University of Michigan)/MED-LauringLab/"
  code_path <- "C:/Users/chbla/Documents/LauringLab/COVID19SequenceCompile/"
  batch_path <- "C:/Users/chbla/Documents/LauringLab/AlertCode"
  
  
} else {
  
  print("User not recognized.")
  
}



################################################################################

options(scipen=999)
source(paste0(code_path, "pipeline_functions.R"))










################################################################################
# SARS-CoV-2 COMPONENT
################################################################################

# sars-cov-2 plate
#<<<<<<< HEAD
plate_name <- "2023_SC2_Nanopore_Run_312"
#=======
plate_name <- "20230605_SC2_Illumina_Run_108"
#>>>>>>> 965c5720bf619c9fbd3a0ac0fa4ba776ee6d272f


################################################################################
#                                 ROUND 1                                      #
################################################################################
# Explanation: Moves proper plate map from load location to use location, compiles
# manifests and plate maps, merges them together, and checks that the plate we're 
# running has a manifest entry for each sample
################################################################################

source(paste0(code_path, "OutsidePipeline/moving_plate_map_files.R"))

source(paste0(code_path, "manifest_code.R"))
# time = 1.06s

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open
# time = 10.56s

source(paste0(code_path, "compile_components_code.R"))
# time = 79.45s

source(paste0(code_path, "OutsidePipeline/subset_compiled_for_fasta.R"))

################################################################################

# checking_sampleids.py step goes here

################################################################################
#                                PANGOLIN                                      #
################################################################################

# Run pangolin in the command line

# Rename this file to plateMapName_pangolin.csv
# Copy this renamed file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/pangolin]
source(paste0(code_path, "OutsidePipeline/moving_pangolin_output.R"))

################################################################################
#                                NEXTCLADE                                     #
################################################################################

# run nextclade from https://clades.nextstrain.org/results

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
# or genbank metadata from ncbi

source(paste0(code_path, "manifest_code.R"))

source(paste0(code_path, "plate_map_code.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "pangolin_code.R"))

source(paste0(code_path, "nextclade_code.R"))

source(paste0(code_path, "gisaid_code.R")) # need to do one final pull, then can retire this
source(paste0(code_path, "genbank_code.R"))

source(paste0(code_path, "compile_components_code.R"))

################################################################################









################################################################################
# INFLUENZA COMPONENT
################################################################################

plate_name <- "20230620_IBV_Illumina_Run_2"
#plate_name <- "DATE_SC2_Illumina_Run_XX_E6440"

plate_name <- "20231107_IAV_Nanopore_Run_40"

num_iavs_on_plate <- 96

################################################################################
### first check, see if any fasta files made for negative control wells
if (grepl("IAV", plate_name)){
  all_fastas <- list.files(paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_name, "/Segment_sequences"))
} else {
  all_fastas <- list.files(paste0(starting_path, "SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/", plate_name, "/Segment_sequences"))
}

if (any(grepl("NC", all_fastas))){
  stop("Negative Control Well generated fasta files.")
}

if (any(grepl("HeLa", all_fastas))){
  stop("HeLa Control Well generated fasta files.")
}

#### if anything prompts a notification, confirm with lab technicians/underworld
#### that plate passed negative control well checks

################################################################################
# ROUND 1
################################################################################

source(paste0(code_path, "OutsidePipeline/moving_plate_map_files.R"))

source(paste0(influenza_path, "influenza_manifestcode.R"))

source(paste0(influenza_path, "influenza_platemapcode.R"))

source(paste0(influenza_path, "influenza_genbankcode.R"))

source(paste0(influenza_path, "influenza_compilecomponentscode.R"))
source(paste0(influenza_path, "influenzaB_compilecomponentscode.R"))

source(paste0(influenza_path, "OutsidePipeline/subset_compiled_influenza.R"))

######
# add in piece here to copy nextclade files
source(paste0(code_path, "OutsidePipeline/moving_nextclade_output.R"))

################################################################################
# ROUND 2
################################################################################

source(paste0(influenza_path, "influenza_manifestcode.R"))

source(paste0(influenza_path, "influenza_platemapcode.R"))

source(paste0(influenza_path, "influenza_nextcladecode.R"))

source(paste0(influenza_path, "influenza_genbankcode.R"))

source(paste0(influenza_path, "influenza_compilecomponentscode.R"))
source(paste0(influenza_path, "influenzaB_compilecomponentscode.R"))

################################################################################
# ROUND 3
################################################################################
source(paste0(influenza_path, "influenza_genbankcode.R"))

source(paste0(influenza_path, "influenza_compilecomponentscode.R"))

################################################################################







################################################################################
# RSVA COMPONENT
################################################################################

plate_name <- "20231005_RSVA_Nanopore_Run_1"
#plate_name <- "20230307_RSVA_Illumina_Run_1_E6444_Nextseq"
#plate_name <- "DATE_SC2_Illumina_Run_XX_E6440"

num_rsvs_on_plate <- 85

################################################################################
# ROUND 1
################################################################################

source(paste0(code_path, "OutsidePipeline/moving_plate_map_files.R"))

source(paste0(code_path, "RSVACode/rsva_manifestcode.R"))

source(paste0(code_path, "RSVACode/rsva_platemapcode.R"))

source(paste0(code_path, "RSVACode/rsva_compilecomponentscode.R"))

source(paste0(code_path, "RSVACode/rsva_subset_compiled_for_fasta.R"))

#                                NEXTCLADE                                     #
################################################################################

# run nextclade from https://clades.nextstrain.org/results

# Copy that re-named file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/nextclade]
source(paste0(code_path, "RSVACode/OutsidePipeline/moving_nextclade_output.R"))

################################################################################
# ROUND 2
################################################################################

source(paste0(code_path, "RSVACode/rsva_manifestcode.R"))

source(paste0(code_path, "RSVACode/rsva_platemapcode.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "RSVACode/rsva_nextclade_code.R"))

source(paste0(code_path, "RSVACode/rsva_compilecomponentscode.R"))

################################################################################
# ROUND 3
################################################################################

source(paste0(code_path, "RSVACode/rsva_compilecomponentscode.R"))








################################################################################
# RSVB COMPONENT
################################################################################

plate_name <- "20231004_RSVB_Nanopore_Run_1"
#plate_name <- "20230307_RSVA_Illumina_Run_1_E6444_Nextseq"
#plate_name <- "DATE_SC2_Illumina_Run_XX_E6440"

num_rsvs_on_plate <- 85

################################################################################
# ROUND 1
################################################################################

source(paste0(code_path, "OutsidePipeline/moving_plate_map_files.R"))


source(paste0(code_path, "RSVBCode/rsvb_manifestcode.R"))

source(paste0(code_path, "RSVBCode/rsvb_platemapcode.R"))

source(paste0(code_path, "RSVBCode/rsvb_compilecomponentscode.R"))


source(paste0(code_path, "RSVBCode/rsvb_subset_compiled_for_fasta.R"))

#                                NEXTCLADE                                     #
################################################################################

# run nextclade from https://clades.nextstrain.org/results

# Copy that re-named file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/nextclade]
source(paste0(code_path, "RSVBCode/OutsidePipeline/moving_nextclade_output.R"))

################################################################################
# ROUND 2
################################################################################

source(paste0(code_path, "RSVBCode/rsvb_manifestcode.R"))

source(paste0(code_path, "RSVBCode/rsvb_platemapcode.R"))
# note: this doesn't run if any plate map excel files are open

source(paste0(code_path, "RSVBCode/rsvb_nextclade_code.R"))

source(paste0(code_path, "RSVBCode/rsvb_compilecomponentscode.R"))

################################################################################
# ROUND 3
################################################################################

source(paste0(code_path, "RSVBCode/rsvb_compilecomponentscode.R"))




