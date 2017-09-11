#!/usr/bin/env Rscript

#The files need to be passed in in order of chromosome

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(ggplot2)
library(plyr)

#Test
#args <- vector(length = 7)
#args[2] <- "/Users/phanstiel4/Documents/code_rep/data/background_20.txt"
#args[3] <- "/Users/phanstiel4/Documents/code_rep/data/TADs_20.txt"
#args[4] <- "/Users/phanstiel4/Documents/code_rep/data/loop_flare_20.txt"
#args[5] <- "/Users/phanstiel4/Documents/code_rep/data/background_3.txt"
#args[6] <- "/Users/phanstiel4/Documents/code_rep/data/TADs_3.txt"
#args[7] <- "/Users/phanstiel4/Documents/code_rep/data/loop_flare_3.txt"


#full data
full_background <- read.delim(args[2], sep = "\t", header = TRUE)

#TADS
tad <- read.delim(args[3], sep = "\t", header = TRUE)

#Loops and flares
landf <- read.delim(args[4], sep = "\t", header = TRUE)

#### Function to read in files ####
pooled_stats <- function(data_df) {
  new_df <- matrix(ncol=5)
  new_df <- as.data.frame(new_df)
  colnames(new_df) <- c("distance", "total", "means","sd", "N")
  
  data_df$sd[is.na(data_df$sd)] <- 0
  data_df$sd <- data_df$sd^2
  for ( i in data_df$distance ) {
    if ( is.na(i) ) next
    if ( i >= 251 ) break
    if ( !(i %in% new_df$distance) ) {
      data_portion <- data_df[data_df$distance == i,]
      data_portion <- na.omit(data_portion)
      if( nrow(data_portion) == 1 ) {
        mm <- data.frame(i, data_portion$total, data_portion$means[1], data_portion$sd[1],data_portion$N[1])
      }
      else {
        total_N <- sum(data_portion$N-1)
        tot <- sum(data_portion$total)
        t <- sum(data_portion$N)
        m <- sum(data_portion$means * ((data_portion$N-1)/total_N))
        v <- sum(data_portion$sd * ((data_portion$N-1)/total_N))
        sd <- sqrt(v)
        mm <- data.frame(i,tot,m,sd,t)
      }
      colnames(mm) <- colnames(new_df)
      new_df <- rbind(new_df, mm)
    }
    
  }
  return(na.omit(new_df))
}

read_in_and_add <- function(file, df) {
  add_df <- read.delim(file, sep = "\t", header = TRUE)
  df <- rbind(df, add_df)
  df <- pooled_stats(df)
  df[is.nan(df$sd)] <- 0
  print(df)
  return(df)
}

#Need to fix pooled stats


### MAIN ###
#Read in the files 3 at a time. The files should be submitted by chromosome in order of full tad and flares
for ( i in seq(5,length(args), by = 3) ) {
  full_background <- read_in_and_add(args[i], full_background)
  tad <- read_in_and_add(args[i+1], tad)
  landf <- read_in_and_add(args[i+2], landf)
}

#Combine
tad$model <- "TADs"
full_background$model <- "Background"
landf$model <- "Loop&FL"

#Diff Mean
tad$tot_mean <- tad$total/tad$N
full_background$tot_mean <- full_background$total/full_background$N
landf$tot_mean <- landf$total/landf$total

combined_data_frame <- rbind(tad,full_background,landf)

p <- ggplot( combined_data_frame[combined_data_frame$distance < 50,], aes( x = distance, y = means, col=model))+
  geom_line()+
  ggtitle("Distance Vs Mean")+
  theme(plot.title = element_text(hjust = .5))

p2 <- ggplot( combined_data_frame[combined_data_frame$distance < 50,], aes( x = distance, y = tot_mean, col=model))+
  geom_line()+
  ggtitle("Distance Vs Mean")+
  theme(plot.title = element_text(hjust = .5))

ggsave(paste(args[1],"/c_dist_vs_sig.png",sep = ""), p)
ggsave(paste(args[1],"/c_dist_vs_sig_tot.png",sep = ""), p2)
write.table(combined_data_frame, file = paste(args[1], "/full_sig_distance_all_things.txt", sep = ""), sep = "\t", row.names = FALSE)
