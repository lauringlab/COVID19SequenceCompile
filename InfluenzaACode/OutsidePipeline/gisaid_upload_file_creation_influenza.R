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

#plate_name <- "20230623_IAV_Illumina_Run_58"

plate_datef <- strsplit(plate_name, "_")[[1]][1] # plate date in YYYYMMDD format
runtech <- strsplit(plate_name, "_")[[1]][3] # nanopore or illumina, will match "PlatePlatform" options
runnum <- strsplit(plate_name, "_")[[1]][5] # number, will match "PlateNumber" options


checking_wd <- getwd()
if (grepl("IAV", plate_name)){
    if (grepl("juliegil", checking_wd)){
      
      #code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
      code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
      starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
      seq_list_path <- paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_datef, "_IAV_", runtech, "_Run_", runnum, "/Segment_sequences/")
      #code_path2 <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/InfluenzaACode/"
      code_path2 <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/InfluenzaACode/"
      
    } else if (grepl("leighbaker", checking_wd)){
      code_path <- "/Users/leighbaker/Documents/Lauring_Lab/COVID19SequenceCompile/"
      starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"
      seq_list_path <- ""
      code_path2 <- ""
      
    } else if (grepl("leighbak", checking_wd)){
      
      starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"
      code_path <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/"
      batch_path <- "/Users/leighbak/Documents/Lauring_Lab/AlertCode"
      seq_list_path <- ""
      code_path2 <- ""
      
    } else {
      
      stop("User not recognized.")
      
    }
} else if (grepl("IBV", plate_name)){
  if (grepl("juliegil", checking_wd)){
    
    #code_path <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/"
    code_path <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/"
    starting_path <- "C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/"
    seq_list_path <- paste0("C:/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/", plate_datef, "_IBV_", runtech, "_Run_", runnum, "/Segment_sequences/")
    #code_path2 <- "/Users/juliegil/Documents/git_synced_code/SequenceCompilationCode/COVID19SequenceCompile/InfluenzaACode/"
    code_path2 <- "C:/Users/juliegil/Documents/UofM_Work/SequenceCompilationCode/InfluenzaACode/"
    
  } else if (grepl("leighbaker", checking_wd)){
    code_path <- "/Users/leighbaker/Documents/Lauring_Lab/COVID19SequenceCompile/"
    starting_path <- "/Users/leighbaker/Dropbox (University of Michigan)/MED-LauringLab/"
    seq_list_path <- ""
    code_path2 <- ""
    
  } else if (grepl("leighbak", checking_wd)){
    
    starting_path <- "/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/"
    code_path <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/"
    batch_path <- "/Users/leighbak/Documents/Lauring_Lab/AlertCode"
    seq_list_path <- ""
    code_path2 <- ""
    
  } else {
    
    stop("User not recognized.")
    
  }
  
  
}

source(paste0(code_path, "pipeline_functions.R"))


################################################################################

if (grepl("IBV", plate_name)){
  dir.create(paste0(starting_path, "/SEQUENCING/INFLUENZA_B/5_GISAID_Uploads/upload_", plate_datef, "_ibv_", tolower(runtech), "_run_", runnum))
  # run comparison code file first, to be sure full_compiled_data matches the one
  # in the secret folder
  source(paste0(code_path2, "OutsidePipeline/comparing_full_secret_influenzaB.R"))
  
  # set output path for gisaid upload file
  # will need to add appropriate folder name at the end of this path
  outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/5_GISAID_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum, "/")
  
  ################################################################################
  
  # read in full compiled pile
  finalfileLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary")
  final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")
  
} else if (grepl("IAV", plate_name)){
  # create gisaid directory
  dir.create(paste0(starting_path, "/SEQUENCING/INFLUENZA_A/5_GISAID_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum))

  ################################################################################
  
  # run comparison code file first, to be sure full_compiled_data matches the one
  # in the secret folder
  source(paste0(code_path2, "OutsidePipeline/comparing_full_secret_influenza.R"))
  
  
  # set output path for gisaid upload file
  # will need to add appropriate folder name at the end of this path
  outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/5_GISAID_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum, "/")
  
  ################################################################################
  
  # read in full compiled pile
  finalfileLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary")
  final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")
  
  
  }

# remove any with missing collection date
final_file <- filter(final_file, !grepl("Replaced with Today Date", flag) & received_source != "WW")

# only keep rows with completeness > 90%
ff <- filter(final_file, as.numeric(nextclade_HA_completeness) >= 90)

# select run of choice
ff <- filter(ff, PlatePlatform == runtech & PlateNumber == runnum)

ff <- filter(ff, subject_id != "")

table(ff$received_source)

################################################################################
# set up alert to duplicate items

if (any(ff$sample_per_subject > 1)){
  print("STOP: Examine this set of GISAID submissions.")
  stop("There are samples from subject_ids that we've sequenced previously.")
}


samples_previous <- filter(ff, sample_per_subject > 1) %>% select(subject_id, sample_id, coll_date)
original_full <- filter(final_file, subject_id %in% unique(samples_previous$subject_id))
### check if the samples are > 90 days apart from one another - then you can let 
### them through.

### uncomment this portion to remove those samples
### to remove these: 
#ff <- filter(ff, sample_per_subject == 1)
#ff <- filter(ff, sample_per_subject == 1 | subject_id == "100432897")
#ff <- filter(ff, sample_id != "007482947")

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
if (any(grepl("IVY", ff$received_source))){
  
  ff$source_state_code <- substr(ff$subject_id, 3, 4)
  cdcivy_manifest_fp <- paste0(starting_path, "SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/Manifests/CDCIVY")
  cdc_sites <- read.csv(paste0(cdcivy_manifest_fp, "/Keys/CDC_SiteCodebook.csv"), colClasses = "character") %>% separate(SiteCode, into = c("site", "state"), sep = "_")
  site_bit <- cdc_sites %>% select(Number, state)
  # add leading zeros to site
  site_bit <- site_bit %>% mutate(Number = case_when(nchar(Number) == 1 ~ paste0("0", Number), 
                                                     T ~ Number))
  #colnames(site_bit) <- c("Number", "state")
  ff <- merge(ff, site_bit, by.x = c("source_state_code"), by.y = c("Number"), all.x = TRUE)
  ff$state <- ifelse(!grepl("IVY", ff$received_source), "", ff$state)
  
} else {
  state <- ""
}

if(any(grepl("RVTN", ff$received_source))){
  
  ff<- ff %>% mutate(state = case_when(grepl("RVTN", received_source) ~ state.abb[match(flag,state.name)], 
                                       T ~ state), 
                     sample_id = case_when(grepl("RVTN", received_source) ~ sample_id_lauring, 
                                           T ~ sample_id))
  
}

if (grepl("IAV", plate_name)){
    ff$StrainName <- ifelse(ff$received_source %in% c("CBR", "UHS", "EPIDIAV", "MFIVE", "HIVE", "MM", "STJ", "UOM", "MARTIN"), paste0("A/Michigan/UOM", ff$sample_id, "/", year(ff$coll_date)), "CHECK")



    ff <- ff %>% mutate(StrainName = case_when(received_source %in% c("CBR", "UHS", "EPIDIAV", "MFIVE", "HIVE", "MM", "STJ", "UOM", "MHOME", "MARTIN") ~ paste0("A/Michigan/UOM", sample_id, "/", year(coll_date)),
                                               received_source %in% c("HFHS", "ASJ") ~ paste0("A/Michigan/MISAPPHIRE", sample_id, "/", year(coll_date)),
                                               grepl("IVY", received_source) ~ paste0("A/", state.name[match(state,state.abb)] , "/IVY", sample_id, "/", year(coll_date)),
                                               grepl("RVTN", received_source) ~ paste0("A/", state.name[match(state,state.abb)] , "/RVTN", sample_id, "/", year(coll_date)),
                                               T ~ "CHECK"
                        ))
    
    ff$Subtype <- ifelse(ff$nextclade_HA_type == "H3", "H3N2", #H3N2
                         ifelse(ff$nextclade_HA_type == "H1", "H1N1", "CHECK")) #H1N1

} else if (grepl("IBV", plate_name)){
  
  ff$StrainName <- ifelse(ff$received_source %in% c("CBR", "UHS", "EPIDIAV", "MFIVE", "HIVE", "MM", "STJ", "UOM", "MARTIN"), paste0("B/Michigan/UOM", ff$sample_id, "/", year(ff$coll_date)), "CHECK")
  
  
  ff <- ff %>% mutate(StrainName = case_when(received_source %in% c("CBR", "UHS", "EPIDIAV", "MFIVE", "HIVE", "MM", "STJ", "UOM", "MHOME", "MARTIN") ~ paste0("B/Michigan/UOM", sample_id, "/", year(coll_date)),
                                             received_source %in% c("HFHS", "ASJ") ~ paste0("B/Michigan/MISAPPHIRE", sample_id, "/", year(coll_date)),
                                             grepl("IVY", received_source) ~ paste0("B/", state.name[match(state,state.abb)] , "/IVY", sample_id, "/", year(coll_date)),
                                             grepl("RVTN", received_source) ~ paste0("B/", state.name[match(state,state.abb)] , "/RVTN", sample_id, "/", year(coll_date)),
                                             T ~ "CHECK"
                                              ))
  
  ff$Subtype <- "B"
  
  
}

if (any(ff$StrainName == "CHECK")){
  stop("Unexpected received source!")
}

#
#ff$Type <- "A"
#ff$SubtypeH <- gsub("H", "", ff$nextclade_HA_type)


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
ff <- ff %>% mutate(OriginatingLabID = case_when(received_source %in% c("CDCIVY", "CDCIVY4", "CDCIVY5", "RVTN") ~ "1960",
                                                 received_source == "HFHS" ~ "3559",
                                                 received_source == "ASJ" ~ "3774",
                                                 received_source == "ASC" ~ "3774",
                                                 T ~ "3201")) # clicical micro lab
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

if (grepl("IAV", plate_name)){
  
    write.csv(ff_gisaid, paste0(starting_path, "SEQUENCING/INFLUENZA_A/5_GISAID_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum, "/gisaid_base.csv"), row.names = FALSE, na = "")

  } else if (grepl("IBV", plate_name)){
    
  write.csv(ff_gisaid, paste0(starting_path, "SEQUENCING/INFLUENZA_B/5_GISAID_Uploads/upload_", plate_datef, "_ibv_", tolower(runtech), "_run_", runnum, "/gisaid_base.csv"), row.names = FALSE, na = "")

  }
#University of Michigan Clinical Microbiology Laboratory
#2800 Plymouth Rd, Ann Arbor, MI, USA



## single upload: A/Michigan/UOM10042526240/2021 (2021-11-21)