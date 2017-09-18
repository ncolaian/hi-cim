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
  dist_vs_counts_tads <- matrix(dist_vs_counts_tads$distance[dist_vs_counts_tads$distance %in% kval], dist_vs_counts_tads$reads[dist_vs_counts_tads$distance %in% kval], ncol = 2)
  flares_and_loops_dvc <- matrix(flares_and_loops_dvc$distance[flares_and_loops_dvc$distance %in% kval], flares_and_loops_dvc$reads[flares_and_loops_dvc$distance %in% kval], ncol = 2)
  background_counts <- matrix(background_counts$distance[background_counts$distance %in% kval], background_counts$reads[background_counts$distance %in% kval], ncol = 2)
  
  #add model and bring df to one thing
  dist_vs_counts_tads$model <- "TADs"
  flares_and_loops_dvc$model <- "Loop&FL"
  background_counts$model <- "Background"
  
  combined <- rbind(dist_vs_counts_tads,flares_and_loops_dvc,background_counts)
  
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
    test3 <- rgamma(length(t$reads), gam_t$estimate[1], rate = gam_t$estimate[2])
    
    
    t_new_ys <- predict(t_spline, seq(1,200))
    t_df <- rbind(t_df, as.data.frame(t_new_ys,col.names = c("distance","means")))
    
    lf_new_ys <- predict(lf_spline, seq(1,200))
    lf_df <- rbind(lf_df, as.data.frame(lf_new_ys,col.names = c("distance","means")))
    
    b_new_ys <- predict(b_spline, seq(1,200))
    b_df <- rbind(b_df, as.data.frame(b_new_ys,col.names = c("distance","means")))
    
    comb2 <- rbind(cbind(t_df, model = "TADs"),cbind(lf_df, model = "Loop&FL"),cbind(b_df, model = "Background" ))
    
    par(mfrow=c(3,3))
    
    g1 <- ggplot( comb_df, aes( x = distance, y = means, col=model))+
      geom_line()+
      ggtitle("Distance Vs Mean - Real")+
      theme(plot.title = element_text(hjust = .5))
    
    g2 <- ggplot( comb2, aes( x = distance, y = means, col=model))+
      geom_line()+
      ggtitle(paste("Smooth Spline Spar =", i, "Distance Vs Mean"))+
      theme(plot.title = element_text(hjust = .5))
    
    print(g1)
    print(g2)
  }
  dev.off()
}