#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)


#Handles input and puts the options passed into the program into a matrix
#Matrix makes it easier to look up values.
params = matrix(c(
  "loop_file", "l", 1, "character",
  "count_file", "cf", 1, "character",
  "chromosome", "chr", 1, "integer",
  "out_dir", "o", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

#test
opt <- c("l", "a", 1)
names(opt) <- c("loop_file","count_file", "chromosome")
opt$loop_file <- "/Users/phanstiel4/Documents/code_rep/data/CI_THP1_O_0.0.0.loops.10Kb.bedpe"
opt$count_file <- "/Users/phanstiel4/Documents/code_rep/data/CI_THP1_O_0.0.0.chr20.10Kb.MQ30.KR_norm.mat"
opt$chromosome <- 20
### Read in data ###

#loop data
loop <- read.delim(opt$loop_file, header = FALSE)
loop_df <- as.data.frame(loop)
colnames(loop_df) <- c("achr", "abin_start", "abin_end", "bchr", "bbin_start", "bbin_end",
                           "bin_length", "loop_pval")
#assign new bin length values to make it a number and not a character
loop_df$distance <- abs(loop_df$abin_start - loop_df$bbin_start)

#norm count data
read_counts <- read.delim(opt$count_file, header = FALSE)
read_counts <- as.data.frame(read_counts)
colnames(read_counts) <- c("start", "end", "reads")

#### Function to create data frames ####
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
    flare_and_loop <- rc_df[rc_df$start >= (start-adjust) & rc_df$end <= (end+adjust),]
    for (j in 1:length(flare_and_loop) ) {
      #get TADS
      if ( (flare_and_loop[j,]$start > (start+adjust)) && (flare_and_loop[j,] < (end-adjust)) ) {
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
  
  #Use ddply to merge all the read counts for a particular distance into a mean. Also report the std
  dist_vs_counts_tads <- ddply(dist_vs_counts_tads, "distance", summarize, means = mean(reads, na.rm = TRUE), sd = sd(reads, na.rm = TRUE), N = sum(!is.na(reads)))
  flares_and_loops_dvc <- ddply(flares_and_loops_dvc, "distance", summarize, means = mean(reads, na.rm = TRUE), sd = sd(reads, na.rm = TRUE), N = sum(!is.na(reads)))
  background_counts <- ddply(background_counts, "distance", summarize, means = mean(reads, na.rm = TRUE), sd = sd(reads, na.rm = TRUE), N = sum(!is.na(reads)))
}

print_out_data <- function(name, dataframe, chr, out) {
  write.table(dataframe, file = paste(out, "/", name, chr, ".txt"), sep = "\t", row.names = FALSE)
}

#### MAIN ####
separate_loops_and_tads(opt$chromosome, loop_df, read_counts)
print_out_data("loop_flare_", flares_and_loops_dvc, opt$chromosome, opt$out_dir)
print_out_data("background_", background_counts, opt$chromosome, opt$out_dir)
print_out_data("TADs_", dist_vs_avg_counts, opt$chromosome, opt$out_dir)

