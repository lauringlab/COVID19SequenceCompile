################################################################################
#                    COVID-19 Sequence Code Functions                          #
#                         Last Updated: 09/01/2022                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

################################################################################
#                   Dates and File Names Formatting                            #
################################################################################

## takes in a file name that has a date as the second "item", separated from the 
## first "item" with a "_"
## Ex. CSTP_YYYYMMDD_1.csv
## checks if the date piece is 8 characters long
## and spits out the date piece in a date format (YYYY-MM-DD)
date_from_file <- function(file_name_in){
  # get the second item out of the file name, as a character string, removing any excess 
  # leading/lagging whitespace
  rec_date <- trimws(as.character(strsplit(file_name_in, "_")[[1]][3]))
  
  # check that the date piece is 8 digits long
  if (nchar(rec_date) != 8){
    print(file_name_in)
    stop("File name date portion (2nd position) is not 8 characters.")
  
  }
  
  # from that piece, pull out the first 4 digits, then the 5th and 6th digits, then the 7th and 8th
  # and place a "-" symbol between those pieces
  rec_date <- paste0(substr(rec_date, 1, 4), "-", substr(rec_date, 5, 6), "-", substr(rec_date, 7, 8))
  
  return(rec_date)
}

## takes in a file name that has a date as the first "item", separated from the 
## other "items" with "_"
## Ex. YYYYMMDD_Plate_#_#_nextclade.tsv
## checks if the date piece is 8 characters long
## and spits out the date piece in a date format (YYYY-MM-DD)
date_from_file_FIRST <- function(file_name_in){
  # get the second item out of the file name, as a character string, removing any excess 
  # leading/lagging whitespace
  rec_date <- trimws(as.character(strsplit(file_name_in, "_")[[1]][1]))
  
  # check that the date piece is 8 digits long
  if (nchar(rec_date) != 8){
    print(file_name_in)
    stop("File name date portion (1st position) is not 8 characters.")
    
  }
  
  # from that piece, pull out the first 4 digits, then the 5th and 6th digits, then the 7th and 8th
  # and place a "-" symbol between those pieces
  rec_date <- paste0(substr(rec_date, 1, 4), "-", substr(rec_date, 5, 6), "-", substr(rec_date, 7, 8))
  
  return(rec_date)
}


### turns the current date into a YYYYMMDD string for use in file name outputs
current_date_string <- function(){
    # add leading zero to month
    if (length(month(Sys.Date())) == 1){
      m <- paste0("0", month(Sys.Date()))
    } else {
      m <- month(Sys.Date())
    }
    # add leading zero to day
    if (length(day(Sys.Date())) == 1){
      d <- paste0("0", day(Sys.Date()))
    } else {
      d <- day(Sys.Date())
    }
    
    today <- paste0(year(Sys.Date()), m, d)
    return(today)
}

################################################################################
#                                Data QA Checks                                #
################################################################################

### takes in a dataframe and source string, and adds a warning to the flag that the 
### subject_id was a "less than" length (assuming from leading zeros)
### also adds in those leading zeros
### input data frame needs to have received_source, subject_id_length columns
subject_id_length_QA <- function(storage_df, source_string){
      
  ### samples from CBR and ED_IDNOW should be Michigan Medicine patients with MRNs
  ### Michigan Medicine MRNs are 9 characters long
  ### samples from CSTP should be UofM affiliated and have UMIDs
  ### UMIDs are 8 characters long
      if (source_string == "CBR" | source_string == "EDIDNOW"){
          # edit flag to note mismatches/instances where leading zeros were re-introduced
          storage_df$flag <- ifelse(is.na(storage_df$flag), "", storage_df$flag)
          storage_df$flag <- ifelse(storage_df$received_source == source_string & storage_df$subject_id_length < 9, 
                                          paste0(storage_df$flag, " ", "MRN < 9 digits + leading 0s restored"), storage_df$flag)
          storage_df$flag <- trimws(storage_df$flag)
          storage_df$flag <- ifelse(storage_df$flag == "", NA, storage_df$flag)
          
          # add in those leading zeros in cases
          
          storage_df$subject_id <- ifelse(storage_df$received_source == source_string & storage_df$subject_id_length < 9, 
                                                with_options(c(scipen = 999), str_pad(storage_df$subject_id, 9, pad = "0")), 
                                                storage_df$subject_id)
      } else if (source_string == "CSTP") {
      
          # edit flag to note mismatches/instances where leading zeros were re-introduced
          storage_df$flag <- ifelse(is.na(storage_df$flag), "", storage_df$flag)
          storage_df$flag <- ifelse(storage_df$received_source == source_string & storage_df$subject_id_length < 8, 
                                          paste0(storage_df$flag, " ", "UMID < 8 digits + leading 0s restored"), storage_df$flag)
          storage_df$flag <- trimws(storage_df$flag)
          storage_df$flag <- ifelse(storage_df$flag == "", NA, storage_df$flag)
          
          # add in those leading zeros in cases
          
          storage_df$subject_id <- ifelse(storage_df$received_source == source_string & storage_df$subject_id_length < 8, 
                                                with_options(c(scipen = 999), str_pad(storage_df$subject_id, 8, pad = "0")), 
                                                storage_df$subject_id)
      } else {
        
          print("Invalid source_string provided")
        
      }
      
      return(storage_df)
}