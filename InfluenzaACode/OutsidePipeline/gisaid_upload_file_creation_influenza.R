################################################################################
#               GISAID File Upload Format Creation for Influenza               #
#                            Created: 11/19/2021                               #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(openxlsx)
library(reshape2)

################################################################################
# just need some of these functions

code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
source(paste0(code_path, "pipeline_functions.R"))

# set starting path
starting_path <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"

################################################################################
### fill in some info manually

plate_datef <- "20220303" # plate date in YYYYMMDD format
runtech <- "Nanopore" # nanopore or illumina, will match "PlatePlatform" options
runnum <- "13" # number, will match "PlateNumber" options

seq_list_path <- paste0("/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_datef, "_IAV_Nanopore_Run_", runnum, "/Segment_sequences/")

################################################################################

# run comparison code file first, to be sure full_compiled_data matches the one
# in the secret folder
code_path2 <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/InfluenzaACode/"
source(paste0(code_path2, "OutsidePipeline/comparing_full_secret_influenza.R"))

# set output path for gisaid upload file
# will need to add appropriate folder name at the end of this path
outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/5_GISAID_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum, "/")

################################################################################

# read in full compiled pile
finalfileLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary")
final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_HA_completeness) >= 90)

# select run of choice
ff <- filter(ff, PlatePlatform == runtech & PlateNumber == runnum)

################################################################################
# set up alert to duplicate items

if (any(ff$sample_per_subject > 1)){
  print("STOP: Examine this set of GISAID submissions.")
  stop("There are samples from subject_ids that we've sequenced previously.")
}


#samples_previous <- filter(ff, sample_per_subject > 1) %>% select(subject_id, sample_id, coll_date)
#original_full <- filter(final_file, subject_id %in% unique(samples_previous$subject_id))
### check if the samples are > 90 days apart from one another - then you can let 
### them through.

### uncomment this portion to remove those samples
### to remove these: 
#ff <- filter(ff, sample_per_subject == 1)
#ff <- filter(ff, sample_per_subject == 1 | sample_id == "10042234896")
#ff <- filter(ff, sample_id != "10041097200")

################################################################################
### fix date formatting
ff <- ff %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

################################################################################

ff$IsolateID <- ""
ff$SegmentIDs <- ""

### strain naming: A/Michigan/UOM[sample_id]/2021
### {flu type A/B} / {collection state} / {UOM}{sample_id} / {collection year}
ff$StrainName <- ifelse(ff$received_source %in% c("CBR", "UHS"), paste0("A/Michigan/UOM", ff$sample_id, "/", year(ff$coll_date)), "CHECK")

if (any(ff$StrainName == "CHECK")){
  stop("Unexpected received source!")
}

#
#ff$Type <- "A"
#ff$SubtypeH <- gsub("H", "", ff$nextclade_HA_type)
ff$Subtype <- ifelse(ff$nextclade_HA_type == "H3", "H3N2", #H3N2
                      ifelse(ff$nextclade_HA_type == "H1", "H1N1", "CHECK")) #H1N1

if(any(ff$Subtype == "CHECK")){
  stop("Issue with subtype assignment")
}

ff$Lineage <- ""
ff$Passage <- "ORI"

ff$Location <- "North America"
ff$Province <- "United States"
ff$SubProvince <- ""
ff$LocationAdditionalInfo <- ""

ff$Host <- "Human"
ff$HostAdditionalInfo <- ""

### now need sequence ids
### read in made file from prep_fasta_gisaid_flu.py

sequence_ids <- read.csv(paste0(seq_list_path, "gisaid_IDlist.csv")) %>% distinct()
sequence_ids <- separate(data = sequence_ids, col = IDS, sep = "\\_", into = c("sample_id", "segment", "platedate"), remove = FALSE)
#sequence_ids <- data.frame(sequence_ids) %>% select(IDS, sample_id, platedate, segment) %>% distinct()
sequence_ids <- reshape2::dcast(sequence_ids, sample_id + platedate ~ segment, value.var = c("IDS"))

sequence_ids <- sequence_ids %>% select(sample_id, HA, `NA`, PB1, PB2, PA, MP, NS, NP)
colnames(sequence_ids) <- c("sample_id", "SeqID_HA", "SeqID_NA", "SeqID_PB1", "SeqID_PB2", "SeqID_PA", "SeqID_MP", "SeqID_NS", "SeqID_NP")

ff <- merge(ff, sequence_ids, by = c("sample_id"))

ff$SeqID_HE <- ""
ff$SeqID_P3 <- ""

ff$SubmittingSampleID <- ""
ff$Authors <- ""
ff$OriginatingLabID <- "3201"
ff$OriginatingSampleID <- ""

ff$CollectionMonth <- ""
ff$CollectionYear <- ""
ff$CollectionDate <- ff$coll_date

ff_gisaid <- ff %>% select(IsolateID, SegmentIDs, StrainName, Subtype, Lineage,
                           Passage, Location, Province, SubProvince, LocationAdditionalInfo, 
                           Host, HostAdditionalInfo, SeqID_HA, SeqID_NA, SeqID_PB1, SeqID_PB2, 
                           SeqID_PA, SeqID_MP, SeqID_NS, SeqID_NP, SeqID_HE, SeqID_P3, 
                           SubmittingSampleID, Authors, OriginatingLabID, OriginatingSampleID, 
                           CollectionMonth, CollectionYear, CollectionDate)

ff_gisaid[is.na(ff_gisaid)] <- ""

write.csv(ff_gisaid, paste0("/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/5_GISAID_Uploads/upload_", plate_datef, "_iav_nanopore_run_", runnum, "/gisaid_base.csv"), row.names = FALSE, na = "")

#University of Michigan Clinical Microbiology Laboratory
#2800 Plymouth Rd, Ann Arbor, MI, USA



## single upload: A/Michigan/UOM10042526240/2021 (2021-11-21)