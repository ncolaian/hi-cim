# This script will be for analyzing the count values vs distance ( Globally and within TADS )

library(plyr)

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
  #Need to implement a test here - Does this look like a tad?
  
}

graph_pot_tads <- ddply(no_test_dist_vs_counts_tads, "distance", summarize, mean = mean(reads), sd = sd(reads))

#plot the signal vs distance curve
plot(dist_vs_avg_counts$distance[dist_vs_avg_counts$distance<50], dist_vs_avg_counts$mean[dist_vs_avg_counts$distance<50],
     xlab = "Distance in 10KB", ylab = "Average Signal", col ="blue", type = "o", main="Signal vs. Distance")
lines(graph_pot_tads$distance[graph_pot_tads$distance < 50], graph_pot_tads$mean[graph_pot_tads$distance < 50], col="red")
legend(20, 1500,c("Global", "Within TADS"), lty = c(1,1), lwd = c(2,2), col=c("blue", "red"))


