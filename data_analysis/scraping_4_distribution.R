#! /usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(plyr)
library(MASS)

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
  dist_vs_counts_tads$used <- sapply(dist_vs_counts_tads, "distance", function(x) {x %in% kval})
  flares_and_loops_dvc$used <- sapply(flares_and_loops_dvc, "distance", function(x) {x %in% kval})
  background_counts$used <- sapply(background_counts, "distance", function(x) {x %in% kval})
    
  #add model and bring df to one thing
  dist_vs_counts_tads$model <- "TADs"
  flares_and_loops_dvc$model <- "Loop&FL"
  background_counts$model <- "Background"
  
  combined <- rbind(dist_vs_counts_tads[dist_vs_counts_tads$used == TRUE,],
                    flares_and_loops_dvc[flares_and_loops$used == TRUE,],
                    background_counts[background_counts$used == TRUE,])
  
  return(combined)
}

print_out_data <- function(name, dataframe, chr, out) {
  write.table(dataframe, file = gsub(" ","",paste(out, "/", name, chr, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE)
}

### MAIN ###
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

#perform stuff
combined_df <- separate_loops_and_tads(opt$chromosome, loop_df, read_counts)
print_out_data("distr", combined_df, opt$chromosome, opt$out_dir)


### This Code is for manual analysis of distributions ###
#***** NEED TO MAKE SURE YOU COMBINE ALL CHROMOSOME DATA ******
fit_distributions <- function(comb_df) {
  dist_vals <- seq(1,60, by=5)
  pdf(file="/Users/phanstiel4/Documents/sim_graphs/distrub_60x5.pdf", w=11, h=8)
  
  for( i in dist_vals) {
    t <- data.frame(comb_df$distance[comb_df$model == "TADs" & comb_df$distance == i], comb_df$reads[comb_df$model == "TADs" & comb_df$distance == i])
    lf <- data.frame(comb_df$distance[comb_df$model == "Loop&FL" & comb_df$distance == i], comb_df$reads[comb_df$model == "Loop&FL" & comb_df$distance == i])
    b <- data.frame(comb_df$distance[comb_df$model == "Background" & comb_df$distance == i], comb_df$reads[comb_df$model == "Background" & comb_df$distance == i])
    colnames(t) <- c("distance", "reads")
    colnames(lf) <- c("distance", "reads")
    colnames(b) <- c("distance", "reads")
    
    #fit the distributions
    nbin_t <- fitdistr(na.omit(as.integer(t$reads)), "negative binomial")
    gam_t <- fitdistr(na.omit(as.integer(t$reads)), "gamma")
    
    nbin_lf <- fitdistr(na.omit(as.integer(t$reads)), "negative binomial")
    gam_lf <- fitdistr(na.omit(as.integer(t$reads)), "gamma")
    
    nbin_b <- fitdistr(na.omit(as.integer(t$reads)), "negative binomial")
    gam_b <- fitdistr(na.omit(as.integer(t$reads)), "gamma")
    
    #Create data from the predicted distributions
    nb_t_df <- rnegbin(length(t$reads), mu = nbin_t$estimate[2], theta = nbin_t$estimate[1])
    gam_t_df <- rgamma(length(t$reads), gam_t$estimate[1], rate = gam_t$estimate[2])
    
    nb_lf_df <- rnegbin(length(lf$reads), mu = nbin_lf$estimate[2], theta = nbin_lf$estimate[1])
    gam_lf_df <- rgamma(length(lf$reads), gam_lf$estimate[1], rate = gam_lf$estimate[2])
    
    nb_b_df <- rnegbin(length(b$reads), mu = nbin_b$estimate[2], theta = nbin_b$estimate[1])
    gam_b_df <- rgamma(length(b$reads), gam_b$estimate[1], rate = gam_b$estimate[2])
    
    t_xlim <- range(0,max(t$distance)+300)
    t_ylim <- range(0, sum(t$reads)/10)
    
    lf_xlim <- range(0,max(t$distance)+300)
    lf_ylim <- range(0, sum(t$reads)/10)
    
    b_xlim <- range(0,max(t$distance)+300)
    b_ylim <- range(0, sum(t$reads)/10)
    
    par(mfcol=c(3,3))
    
    t1 <- hist(t$reads, breaks = 200, xlim = t_xlim, ylim = t_ylim, main = paste("Original Read Distribution Dist =", i) )
    t2 <- hist(nb_t_df, breaks = 200, xlim = t_xlim, ylim = t_ylim, main = paste("Simulated Neg.Binomial Distribution =", i))
    t3 <- hist(gam_t_df, breaks = 200, xlim = t_xlim, ylim = t_ylim, main = paste("Simulated Gamma Distribution =", i))
    
    lf1 <- hist(lf$reads, breaks = 200, xlim = lf_xlim, ylim = lf_ylim, main = paste("Original Read Distribution =", i))
    lf2 <- hist(nb_lf_df, breaks = 200, xlim = lf_xlim, ylim = lf_ylim, main = paste("Simulated Neg.Binomial Distribution =", i))
    lf3 <- hist(gam_lf_df, breaks = 200, xlim = lf_xlim, ylim = lf_ylim, main = paste("Simulated Gamma Distribution =", i))
    
    b1 <- hist(b$reads, breaks = 200, xlim = b_xlim, ylim = b_ylim, main = paste("Original Read Distribution =", i))
    b2 <- hist(nb_b_df, breaks = 200, xlim = b_xlim, ylim = b_ylim, main = paste("Simulated Neg.Binomial Distribution =", i))
    b3 <- hist(gam_b_df, breaks = 200, xlim = b_xlim, ylim = b_ylim, main = paste("Simulated Gamma Distribution =", i))
    
    print(t1,t2,t3)
    print(lf1,lf2,lf3)
    print(b1,b2,b3)
  }
  dev.off()
}