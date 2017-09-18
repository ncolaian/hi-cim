#! /usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(plyr)

params = matrix(c(
  "loop_file", "l", 1, "character",
  "count_file", "cf", 1, "character",
  "chromosome", "chr", 1, "integer",
  "out_dir", "o", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)


separate_loops_and_tads <- function(chr, loops, rc_df) {
  #Set up tad data frame
  dist_vs_counts_tads <- matrix(ncol = 2)
  dist_vs_counts_tads <-as.data.frame(dist_vs_counts_tads)
  colnames(dist_vs_counts_tads) <- c("distance", "reads")
  
  #Set up loop data frame
  flares_and_loops_dvc <- matrix(ncol=2)
  flares_and_loops_dvc <- as.data.frame(flares_and_loops_dvc)
  colnames(flares_and_loops_dvc) <- c("distance", "reads")
  
  #set up vector to keep track of locations in a loop or flare
  data_used_vec <- c()
  chr <- paste("chr",chr, sep = "")
  #go through each loop from the chromosome given
  for ( i in 1:length(loop_df$achr[loop_df$achr == chr])) {
    #get individual loops
    lp_info <- loop_df[loop_df$achr == chr,][i,]
    #need to get all the reads from bins equal to this and within it
    start <- (lp_info$abin_start/10000)
    end <- (lp_info$bbin_start/10000)
    adjust <- 2
    #This gets a big tad containinf tad and flares
    flare_and_loop <- rc_df[rc_df$start >= (start-adjust) & rc_df$end <= (end+adjust),]
    for (j in 1:length(flare_and_loop) ) {
      if ( (is.na(flare_and_loop[j,]$start))) next;
      #get TADS  ***Edited this recently***
      if ( (flare_and_loop[j,]$start > (start+adjust)) && (flare_and_loop[j,]$end < (end-adjust)) ) {
        mm_matrix <- data.frame((flare_and_loop[j,]$end - flare_and_loop[j,]$start), 
                                flare_and_loop[j,]$reads)
        colnames(mm_matrix) <- colnames(dist_vs_counts_tads)
        dist_vs_counts_tads <- rbind(dist_vs_counts_tads, mm_matrix)
      }
      #get Flares
      else {
        mm_matrix <- data.frame((flare_and_loop[j,]$end - flare_and_loop[j,]$start), 
                                flare_and_loop[j,]$reads)
        colnames(mm_matrix) <- colnames(flares_and_loops_dvc)
        flares_and_loops_dvc <- rbind(flares_and_loops_dvc, mm_matrix)
      }
      data_used_vec <- c(data_used_vec, (paste(as.character(flare_and_loop[j,]$start), ",", as.character(flare_and_loop[j,]$end))))
    }
  }
  #Create a data frame with all the data not in a tad or flare/loop
  data_used_vec <- unique(data_used_vec)
  background_counts <- rc_df[!(paste(as.character(rc_df$start), ",", as.character(rc_df$end)) %in% data_used_vec),]
  background_counts <- data.frame(abs(background_counts$end - background_counts$start),background_counts$reads)
  colnames(background_counts) <- c("distance", "reads")
  
  #only keep the the reads that are in spots of 5
  kval <- seq(1,60,by = 5)
  dist_vs_counts_tads <- matrix(dist_vs_counts_tads$distance[dist_vs_counts_tads$distance %in% kval], dist_vs_counts_tads$reads[dist_vs_counts_tads$distance %in% kval], ncol = 2)
  flares_and_loops_dvc <- matrix(flares_and_loops_dvc$distance[flares_and_loops_dvc$distance %in% kval], flares_and_loops_dvc$reads[flares_and_loops_dvc$distance %in% kval], ncol = 2)
  background_counts <- matrix(background_counts$distance[background_counts$distance %in% kval], background_counts$reads[background_counts$distance %in% kval], ncol = 2)
  
  return(list(flares_and_loops_dvc, background_counts, dist_vs_counts_tads))
}

print_out_data <- function(name, dataframe, chr, out) {
  write.table(dataframe, file = gsub(" ","",paste(out, "/", name, chr, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE)
}

### MAIN ###
data_frame_list <- separate_loops_and_tads(opt$chromosome, loop_df, read_counts)

print_out_data("loop_flare_", data_frame_list[[1]], opt$chromosome, opt$out_dir)
print_out_data("background_", data_frame_list[[2]], opt$chromosome, opt$out_dir)
print_out_data("TADs_", data_frame_list[[3]], opt$chromosome, opt$out_dir)