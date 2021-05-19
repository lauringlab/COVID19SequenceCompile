################################################################################
#                    COVID-19 Sequence Code Functions                          #
#                         Last Updated: 05/19/2021                             #
#                 Code Edited By: Julie (Jules) Gilbert                        #
################################################################################

## takes in a file name that has a date as the second "item", separated from the 
## first "item" with a "_"
## Ex. CSTP_YYYYMMDD_1.csv
## checks if the date piece is 8 characters long
## and spits out the date piece in a date format (YYYY-MM-DD)
date_from_file <- function(file_name_in){
  # get the second item out of the file name, as a character string, removing any excess 
  # leading/lagging whitespace
  rec_date <- trimws(as.character(strsplit(file_name_in, "_")[[1]][2]))
  
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
