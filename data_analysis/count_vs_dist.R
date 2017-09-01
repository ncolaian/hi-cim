# This script will b be for analyzing the count values vs distance ( Globally and within TADS )

library(plyr)
library(ggplot2)

### Read in data ###

#loop data
loop <- read.delim("/Users/phanstiel4/Documents/code_rep/data/CI_THP1_O_0.0.0.loops.10Kb.bedpe", header = FALSE)
loop_matrix <- as.data.frame(loop)
colnames(loop_matrix) <- c("achr", "abin_start", "abin_end", "bchr", "bbin_start", "bbin_end",
                           "bin_length", "loop_pval")
#norm count data
read_counts <- read.delim("/Users/phanstiel4/Documents/code_rep/data/CI_THP1_O_0.0.0.chr20.10Kb.MQ30.KR_norm.mat", header = FALSE)
read_counts <- as.data.frame(read_counts)
colnames(read_counts) <- c("start", "end", "reads")

#assign new bin length values to make it a number and not a character
loop_matrix$distance <- abs(loop_matrix$abin_start - loop_matrix$bbin_start)

#create a distance to count matrix
dist_vs_avg_counts <- data.frame(abs(read_counts$start-read_counts$end), read_counts$reads)
colnames(dist_vs_avg_counts) <- c("distance", "reads")

dist_vs_avg_counts <- ddply(dist_vs_avg_counts, "distance", summarize, mean = mean(reads), sd = sd(reads))

#find domains and find the signal vs distance ratio in domains with loops
loop_matrix$achr <- as.character(loop_matrix$achr)
loop_matrix$bchr <- as.character(loop_matrix$bchr)

#### GLOBAL DISTANCE VS Signal CODE ###
no_test_dist_vs_counts_tads <- matrix(ncol = 2)
no_test_dist_vs_counts_tads <-as.data.frame(no_test_dist_vs_counts_tads)
colnames(no_test_dist_vs_counts_tads) <- c("distance", "reads")
for ( i in 1:length(loop_matrix$achr[loop_matrix$achr == "chr20"])) {
  #get individual loops
  lp_info <- loop_matrix[loop_matrix$achr == "chr20",][i,]
  
  #need to get all the reads from bins equal to this and within it
  start <- (lp_info$abin_start/10000)
  end <- (lp_info$bbin_start/10000)
  
  for (j in 1:length(read_counts[read_counts$start >= start & read_counts$end <= end,]) ) {
    mm_matrix <- data.frame((read_counts[read_counts$start > start & read_counts$end < end,][j,]$end - 
                        read_counts[read_counts$start > start & read_counts$end < end,][j,]$start), 
                        read_counts[read_counts$start > start & read_counts$end < end,][j,]$reads)
    colnames(mm_matrix) <- colnames(no_test_dist_vs_counts_tads)
    no_test_dist_vs_counts_tads <- rbind(no_test_dist_vs_counts_tads, mm_matrix)
  }
}

graph_pot_tads <- ddply(no_test_dist_vs_counts_tads, "distance", summarize, means = mean(reads), sd = sd(reads))

#plot the signal vs distance curve
plot(dist_vs_avg_counts$distance[dist_vs_avg_counts$distance<50], dist_vs_avg_counts$mean[dist_vs_avg_counts$distance<50],
     xlab = "Distance in 10KB", ylab = "Average Signal", col ="blue", type = "o", main="Signal vs. Distance")
lines(graph_pot_tads$distance[graph_pot_tads$distance < 50], graph_pot_tads$mean[graph_pot_tads$distance < 50], col="red")
legend(20, 1500,c("Global", "Within TADS"), lty = c(1,1), lwd = c(2,2), col=c("blue", "red"))

#Need to implement a test here - Does this look like a tad?
#using an implementation of the arrowhead method that the Rao paper described (Thanks to Katie I sorta understand it)
######### TAD TESTING CODE ############
for ( i in 1:length(loop_matrix$achr[loop_matrix$achr == "chr20"])) {
  #get individual loops
  lp_info <- loop_matrix[loop_matrix$achr == "chr20",][i,]
  
  #need to get all the reads from bins equal to this and within it
  start <- (lp_info$abin_start/10000)
  end <- (lp_info$bbin_start/10000)
  
  #this returns a list of vectors for the outer(first two) and inner(last 2) positions of comparison
  bot_and_top_vecs <- arrowhead_location_vecs(start, end)
  
  #create an outer and inner value and perform a t-test on this.

}


## Create Location Vectors ##
#This creates a list of 4 vectors that represent the locations where we want to consider read counts
arrowhead_location_vecs <- function( starting, ending ) {
  d <- ending - starting
  half_d <- trunc(d/2)
  #create outer top and outer bot
  top_max_h <- ending + d
  bot_min_w <- starting - d
  out_bot_vec <- c()
  out_top_vec <- c()
  
  marker <- 0
  add <- 1
  for ( i in 1:(d-1) ) {
    #create a stepping system
    if ( marker == 2 ) {
      add <- add+1
      marker <- 0
    }
    marker <- marker + 1
    
    #for the outs
    for ( j in 1:add ) {
      out_bot_vec <- c(out_bot_vec, (paste(as.character(bot_min_w+i), ",", as.character(starting+j))))
      out_top_vec <- c(out_top_vec, (paste(as.character(ending-j), ",", as.character(top_max_h-i)))) 
    }
  }
  
  #create the inner top and inner bot 
  in_top_vec <- c()
  in_bot_vec <- c()
  for (i in 1:length(out_top_vec)) {
    num_string_top <- out_top_vec[i]
    num_string_bot <- out_bot_vec[i]
    
    top_vec <- strsplit(num_string_top, " , ")
    bot_vec <- strsplit(num_string_bot, " , ")
    
    d <- as.integer(top_vec[[1]][2]) - as.integer(top_vec[[1]][1])  #distances are the same for each
    
    in_top_vec <- c(in_top_vec, (paste(as.character(as.integer(top_vec[[1]][1]) - d), ",", top_vec[[1]][1])))
    in_bot_vec <- c(in_bot_vec, (paste(bot_vec[[1]][2], ",", as.character(as.integer(bot_vec[[1]][2])+d))))
  }
  
  return(list(out_top_vec, out_bot_vec, in_top_vec, in_bot_vec) )
}

#test
arrowhead_location_vecs(10,15)

#linear regression model
tail(dist_vs_avg_counts)
m<- lm(mean ~ distance, data = dist_vs_avg_counts[dist_vs_avg_counts$distance<200,])
summary(m)

m <- nls( mean ~ a/(distance+1), data = dist_vs_avg_counts[dist_vs_avg_counts$distance<20,], start=list(a=2000))
summary(m)

ggplot( dist_vs_avg_counts[dist_vs_avg_counts$distance<20,], aes( x = distance, y= mean ))+
  geom_point()+
  geom_smooth(method = "nls", formula = y ~ a/(x+1),se = F, method.args = list(start = c(a = 3000)))



