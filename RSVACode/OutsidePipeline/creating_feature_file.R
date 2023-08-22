# rsv a feature table creation

location1 <- "/Users/juliegil/Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/RSV_A/6_GenBank_Uploads/FEATURE_TABLES/test"

full_list_in <- read.csv(paste0(location1, "/RSVA_codon_positions.csv"))

### need to melt this dataset

full_melt <- reshape2::melt(full_list_in, id.vars = c("gene", "position"))

# replace . in sample names with "-"
full_melt <- full_melt %>% mutate(variable = gsub("\\.", "-", variable))


set_of_lines <- c()
### for every sample in full_melt, i need to make a txt file formatted entry
### all concatenated into one set
for (every_sample in unique(full_melt$variable)){
    one_sample <- filter(full_melt, variable == every_sample) %>% arrange(value)
    
    # first write the header 
    
    set_of_lines <- c(set_of_lines, paste0(">", unique(one_sample$variable)))
    
    
    # move down rows
    for (i in seq(1, nrow(one_sample))){
      # gather a line
      if ((i %% 2) == 0){
        # even
        # then grab end value and add CDS tab and Product space segment
        aline <- paste0(aline, one_sample[i, 4], "\tCDS\tProduct ", one_sample[i, 1])
        # then write that line
        set_of_lines <- c(set_of_lines, aline)
      } else {
        # odd
        aline <- "" # set new line
        # only need to grab value
        aline <- paste0(aline, one_sample[i, 4], "\t")
      }
      
      
    }
}

fileConn <- file(paste0(location1, "/output.txt"))
writeLines(set_of_lines, fileConn)
close(fileConn)
