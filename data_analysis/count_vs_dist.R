# This script will be for analyzing the count values vs distance ( Globally and within TADS )

library(plyr)

### Read in data ###

#loop data
loop_matrix <- read.delim(, header = FALSE)
loop_matrix <- as.data.frame(loop_matrix)
colnames(loop_matrix) <- c("1chr", "1bin_start", "1bin_end", "2chr", "2bin_start", "2bin_end",
                           "bin_length", "loop_pval")
#norm count data
read_counts <- read.delim(, header = FALSE)
read_counts <- as.data.frame(read_counts)
colnames(read_counts) <- c("start", "end", "reads")

#assign new bin length values to make it a number and not a character
loop_matrix$bin_length <- abs(loop_matrix$1bin_start - loop_matrix$2bin_start)

#create a distance to count matrix
dist_vs_avg_counts <- data.frame(c((read_counts$start-read_counts$end), read_counts$reads) ncol=2)
colnames(dist_vs_avg_counts) <- c("distance", "reads")

dist_vs_avg_counts <- ddply(dist_vs_avg_counts, "reads", summarize, mean = mean(reads))

