################################################################################
#               GenBank File Upload Format Creation for Influenza              #
#                            Created: 9/6/2023                                 #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

library(tidyverse)
library(openxlsx)
library(reshape2)

################################################################################
# just need some of these functions

#plate_name <- "20240130_IBV_Nanopore_Run_5"

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
    seq_list_path <- paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_datef, "_IAV_", runtech, "_Run_", runnum, "/Segment_sequences/")
    code_path2 <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/InfluenzaACode/"
    
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
    seq_list_path <-paste0("/Users/leighbak/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/", plate_datef, "_IBV_", runtech, "_Run_", runnum, "/Segment_sequences/")
    code_path2 <- "/Users/leighbak/Documents/Lauring_Lab/COVID19SequenceCompile/InfluenzaACode/"
    
  } else {
    
    stop("User not recognized.")
    
  }
  
  
}

source(paste0(code_path, "pipeline_functions.R"))

#print(starting_path)
################################################################################


################################################################################

if (grepl("IBV", plate_name)){
  dir.create(paste0(starting_path, "/SEQUENCING/INFLUENZA_B/6_GenBank_uploads/upload_", plate_datef, "_ibv_", tolower(runtech), "_run_", runnum))
  # run comparison code file first, to be sure full_compiled_data matches the one
  # in the secret folder
  source(paste0(code_path2, "OutsidePipeline/comparing_full_secret_influenzaB.R"))
  
  # set output path for gisaid upload file
  # will need to add appropriate folder name at the end of this path
  outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/6_GenBank_uploads/upload_", plate_datef, "_ibv_", tolower(runtech), "_run_", runnum, "/")
  
  ################################################################################
  
  # read in full compiled pile
  finalfileLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary")
  final_file <- read.csv(paste0(finalfileLOC, "/full_compiled_data.csv"), colClasses = "character")
  
} else if (grepl("IAV", plate_name)){
  # create gisaid directory
  dir.create(paste0(starting_path, "/SEQUENCING/INFLUENZA_A/8_GenBank_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum))
  
  ################################################################################
  
  # run comparison code file first, to be sure full_compiled_data matches the one
  # in the secret folder
  source(paste0(code_path2, "OutsidePipeline/comparing_full_secret_influenza.R"))
  
  
  # set output path for gisaid upload file
  # will need to add appropriate folder name at the end of this path
  outputLOC <- paste0(starting_path, "SEQUENCING/INFLUENZA_A/8_GenBank_Uploads/upload_", plate_datef, "_iav_", tolower(runtech), "_run_", runnum, "/")
  
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
  print("STOP: Examine this set of Genbank submissions.")
  stop("There are samples from subject_ids that we've sequenced previously.")
}


samples_previous <- filter(ff, sample_per_subject > 1) %>% select(subject_id, sample_id, coll_date)
original_full <- filter(final_file, subject_id %in% unique(samples_previous$subject_id))
### check if the samples are > 90 days apart from one another - then you can let 
### them through.

### uncomment this portion to remove those samples
### to remove these: 
#ff <- filter(ff, sample_per_subject == 1)
#ff <- filter(ff, sample_per_subject == 1 | subject_id %in% c("045688124", "101522622"))
#ff <- filter(ff, sample_id != "007482947")

################################################################################
### fix date formatting
ff <- ff %>% mutate(coll_date = case_when(grepl("/", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%m/%d/%Y")), 
                                          grepl("-", coll_date) ~ as.character(as.POSIXct(coll_date, format = "%Y-%m-%d")), 
                                          T ~ NA_character_))

################################################################################

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
                                       T ~ state)#, 
                     #sample_id = case_when(grepl("RVTN", received_source) ~ sample_id_lauring, 
                                           #T ~ sample_id)
                     )
  
}


#ff$State <- state.name[match(ff$StateAbbrev,state.abb)]
ff <- ff %>% mutate(StateAbbrev = case_when(received_source == "RVTN" ~ state, 
                                      received_source == "CDCIVY" ~ state, 
                                      T ~ "MI"))

ff <- ff %>% mutate(State = state.name[match(StateAbbrev,state.abb)])

ff <- ff %>% mutate(Location = paste0("USA: ", State))



# creating VirusName which will be the Sequence_ID

ff <- ff %>% mutate(VirusName = case_when(received_source == "CDCIVY" ~ paste0("IVY-", sample_id),
                                          received_source == "CDCIVY4" ~ paste0("IVY-", sample_id),
                                          received_source == "CDCIVY5" ~ paste0("IVY-", sample_id),
                                          received_source == "CDCIVY6" ~ paste0("IVY-", sample_id),
                                          received_source == "RVTN" ~ paste0("RVTN-", sample_id_lauring),
                                          received_source == "VIEW" ~ paste0("VIEW-", sample_id_lauring),
                                          received_source == "IVYIC" ~ paste0("IVYIC-", sample_id),
                                          received_source == "MDHHS" ~ paste0("UM-", sample_id),
                                          received_source == "CBR" ~ paste0("UM-", sample_id),
                                          received_source == "HFHS" ~ paste0("MIS-", sample_id),
                                          received_source == "ASC" ~ paste0("MIS-", sample_id),
                                          received_source == "ASJ" ~ paste0("MIS-", sample_id),
                                          received_source == "TRINITY" ~ paste0("MIS-", sample_id),
                                          T ~ paste0("UM-", sample_id)))


 ##create genbank formatted isolate name
 # ff <- ff %>% mutate(IsolateName = case_when(received_source == "CDCIVY" ~ paste0("SARS-CoV-2/human/USA/", "IVY-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "CDCIVY4" ~ paste0("SARS-CoV-2/human/USA/", "IVY-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "CDCIVY5" ~ paste0("SARS-CoV-2/human/USA/", "IVY-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "CDCIVY6" ~ paste0("SARS-CoV-2/human/USA/", "IVY-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "RVTN" ~ paste0("SARS-CoV-2/human/USA/", "RVTN-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "VIEW" ~ paste0("SARS-CoV-2/human/USA/", "VIEW-", sample_id_lauring, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "IVYIC" ~ paste0("SARS-CoV-2/human/USA/", "IVYIC-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "HFHS" ~ paste0("SARS-CoV-2/human/USA/", "MIS-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "ASC" ~ paste0("SARS-CoV-2/human/USA/", "MIS-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "ASJ" ~ paste0("SARS-CoV-2/human/USA/", "MIS-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "TRINITY" ~ paste0("SARS-CoV-2/human/USA/", "MIS-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           received_source == "MDHHS" ~ paste0("SARS-CoV-2/human/USA/", "MI-UM-", sample_id, "/", substr(coll_date, 1, 4)),
 #                                           T ~ paste0("SARS-CoV-2/human/USA/MI-UM-", sample_id, "/", substr(coll_date, 1, 4))))
 # #received_source == "CDCIVY4" ~ paste0("SARS-CoV-2/human/USA/", StateAbbrev, "/", "IVY-", sample_id, "/", substr(coll_date, 1, 4)),

if (any(nchar(ff$VirusName) > 24)){
  print(filter(ff, VirusName > 24))
  stop("Virus name has too many characters for GenBank.")
}


############### filter down to only things that have 8 segments
all_files <- list.files(seq_list_path, pattern = "*.fasta")

check_segment_number <- data.frame(all_files)
check_segment_number <- check_segment_number %>% separate(all_files, c('first', 'second'), sep = "\\.", remove = FALSE)
check_segment_number <- check_segment_number %>% separate(first, c('sample_id', 'type', 'segment1', 'segment2'), sep = "_", remove = FALSE)
check_segment_number <- check_segment_number %>% group_by(sample_id) %>% mutate(count = length(unique(segment1)))
check_segment_number <- filter(check_segment_number, count == 8)

ff <- merge(ff, check_segment_number, by = c("sample_id"))
ff[is.na(ff)] <- ""

if (grepl("IAV", plate_name)){

ff <- ff %>% mutate(Sequence_ID = trimws(paste0(VirusName, segment1, segment2)), # needs to be virus name + sequence name ids 
                    organism = "Influenza A virus",
                    isolate = VirusName, 
                    country = paste0("USA:", State), 
                    host = "Homo sapiens", 
                    collection_date = as.character(coll_date), 
                    isolation_source = "patient isolate", 
                    serotype = case_when(nextclade_HA_type == "H1" ~ "H1N1", 
                                         nextclade_HA_type == "H3" ~ "H3N2", 
                                         T ~ "unknown"))

} else {
  ff <- ff %>% mutate(Sequence_ID = trimws(paste0(VirusName, segment1, segment2)), # needs to be virus name + sequence name ids 
                      organism = "Influenza B virus",
                      isolate = VirusName, 
                      country = paste0("USA:", State), 
                      host = "Homo sapiens", 
                      collection_date = as.character(coll_date), 
                      isolation_source = "patient isolate", 
                      serotype = case_when(nextclade_HA_type == "victoria" ~ "Victoria", 
                                           nextclade_HA_type == "yamagata" ~ "Yamagata", 
                                           T ~ "unknown"))
  
}

if (any(ff$serotype == "unknown")){
  stop("Unknown serotype")
}

# Oxford Nanopore, Illumina MiSeq
# ff$SequencingTechnology <- ifelse(ff$PlatePlatform == "Nanopore", "Oxford Nanopore Midnight", 
#                                   ifelse(ff$PlatePlatform == "Illumina", "Illumina NextSeq 1000", "Unknown"))
# 
# unknown_tech <- filter(ff, SequencingTechnology == "Unknown")
# 
# if (nrow(unknown_tech) != 0){
#   stop("Check Sequencing Technology options.")
# }
# 
# ### Assembly Method
# ff$AssemblyMethod <- ifelse(ff$PlatePlatform == "Nanopore", "ARTIC Network pipeline Midnight", 
#                             ifelse(ff$PlatePlatform == "Illumina", "BWA-MEM, iVar", "Unknown"))
# 
# unknown_assembly <- filter(ff, AssemblyMethod == "Unknown")

# if (nrow(unknown_assembly) != 0){
#   stop("Check Assembly Method options.")
# }
# 

################################################################################

### write out VirusName + sample_id crosswalk for use in making 
# .all.consensus.final.genbank.fasta

ff_crosswalk <- ff %>% select(sample_id, VirusName) %>% distinct()

if (grepl("IAV", plate_name)){

    write.csv(ff_crosswalk, paste0(starting_path, "SEQUENCING/INFLUENZA_A/3_ProcessedGenomes/", plate_datef, "_IAV_", runtech, "_Run_", runnum, "/", plate_datef, "_IAV_", runtech, "_Run_", runnum, ".forgenbank.meta.csv"), row.names = FALSE, na = "")

} else {
  
  write.csv(ff_crosswalk, paste0(starting_path, "SEQUENCING/INFLUENZA_B/3_ProcessedGenomes/", plate_datef, "_IBV_", runtech, "_Run_", runnum, "/", plate_datef, "_IBV_", runtech, "_Run_", runnum, ".forgenbank.meta.csv"), row.names = FALSE, na = "")  

}

if (grepl("IBV", plate_name)){
## select variables
ff_writeout <- ff %>% select(Sequence_ID, organism, country, host, isolate, collection_date, isolation_source)#, serotype)
colnames(ff_writeout) <- c("Sequence_ID", "Organism", "country", "host", "isolate", "collection-date", "isolation-source")#, "serotype")
} else {
  ff_writeout <- ff %>% select(Sequence_ID, organism, country, host, isolate, collection_date, isolation_source, serotype)
  colnames(ff_writeout) <- c("Sequence_ID", "Organism", "country", "host", "isolate", "collection-date", "isolation-source", "serotype")
  
}
# Sequence_ID	Organism	country	host	isolate	collection-date	isolation-source	serotype
ff_writeout <- ff_writeout %>% distinct()

## genbank upload file name
today <- current_date_string()
gufn <- paste0(today, "_Lauring_genbank_upload_metadata_run_", runnum) 


write.table(ff_writeout, paste0(outputLOC, gufn, ".tsv"), sep = "\t", row.names = FALSE, na = "", )

# ## write to excel file (follow format)
# wb <- loadWorkbook(paste0(starting_path, "/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/GISAID_UPLOAD_TEMPLATE2.xlsx"))
# 
# # fill in the submissions tab with built data frame
# writeData(wb, ff_writeout, sheet = "Submissions", startRow = 3, startCol = 1, colNames = FALSE)
# 
# saveWorkbook(wb, paste0(outputLOC, gufn, ".xlsx"), overwrite = TRUE)

