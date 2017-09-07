#!/usr/bin/env Rscript

#The files need to be passed in in order of chromosome

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(ggplot2)
library(plyr)

#full data
full_background <- matrix(ncol=4)
full_background <- as.data.frame(full_background)
colnames(full_background) <- c("distance", "means", "sd", "N")

#TADS
tad <- matrix(ncol=4)
tad <- as.data.frame(tad)
colnames(tad) <- c("distance", "means", "sd", "N")

#Loops and flares
landf <- matrix(ncol=4)
landf <- as.data.frame(landf)
colnames(landf) <- c("distance", "means","sd", "N")

#### Function to read in files ####
pooled_stats <- function(data) {
  new_df <- matrix(col=4)
  new_df <- as.data.frame(new_df)
  colnames(new_df) <- c("distance", "means","sd", "N")
  
  data$sd <- data$sd^2
  for ( i in data$distance ) {
    if ( !(i %in% new_df$distance) ) {
      distance_portion <- data[i %in% data$distance,]
      total_N <- sum(data_portion$N-1)
      t <- sum(data_portion$N)
      m <- sum(data_portion$means * ((data_portion$N-1)/total_N))
      v <- sum(data_portion$sd * ((data_portion$N-2)/total_N))
      sd <- sqrt(v)
      mm <- data.frame(i,m,sd,t)
      colnames(mm) <- colnames(new_df)
      rbind(new_df, mm)
    }
  }
  
  return(new_df)
}

read_in_and_add <- function(file, df) {
  add_df <- read.delim(file, sep = "\t", header = TRUE)
  df <- rbind(df, add_df)
}




### MAIN ###
#Read in the files 3 at a time. The files should be submitted by chromosome in order of full tad and flares
for ( i in seq(1,3,length(args)) ) {
  read_in_and_add(args[i], full_background)
  read_in_and_add(args[i+1], tad)
  read_in_and_add(args[i+2], landf)
}

#create pooled stats
full_background <- pooled_stats(full_background)
tad <- pooled_stats(tad)
landf <- pooled_stats(landf)

#Combine
tad$model <- "TADs"
full_background$model <- "Background"
landf$model <- "Loop&FL"

combined_data_frame <- rbind(tad,full_background,landf)

p <- ggplot( combined_data_frame[combined_data_frame$distance < 50,], aes( x = distance, y = means, col=model))+
  geom_line()+
  ggtitle("Distance Vs Mean")+
  theme(plot.title = element_text(hjust = .5))

ggsave("c_dist_vs_sig", p)
