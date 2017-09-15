# This script will b be for analyzing the count values vs distance ( Globally and within TADS )

library(plyr)
library(ggplot2)
library(vcd)
library(MASS)

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

flares_and_loops_dvc <- matrix(ncol=2)
flares_and_loops_dvc <- as.data.frame(flares_and_loops_dvc)
colnames(flares_and_loops_dvc) <- c("distance", "reads")


for ( i in 1:length(loop_matrix$achr[loop_matrix$achr == "chr20"])) {
  #get individual loops
  lp_info <- loop_matrix[loop_matrix$achr == "chr20",][i,]
  
  #need to get all the reads from bins equal to this and within it
  start <- (lp_info$abin_start/10000)
  end <- (lp_info$bbin_start/10000)
  adjust <- 2
  data_used_vec <- c()
  flare_and_loop <- read_counts[read_counts$start >= (start-adjust) & read_counts$end <= (end+adjust),]
  for (j in 1:length(flare_and_loop) ) {
    #get TADS
    if ( flare_and_loop[j,]$start > (start+adjust) && flare_and_loop[j,]$end < (end-adjust) ) {
      mm_matrix <- data.frame((flare_and_loop[j,]$end - flare_and_loop[j,]$start), 
                              flare_and_loop[j,]$reads)
      colnames(mm_matrix) <- colnames(no_test_dist_vs_counts_tads)
      no_test_dist_vs_counts_tads <- rbind(no_test_dist_vs_counts_tads, mm_matrix)
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
trial_counts <- read_counts[!(paste(as.character(read_counts$start), ",", as.character(read_counts$end)) %in% data_used_vec),]
no_tad_or_loops <- data.frame(abs(trial_counts$end-trial_counts$start),trial_counts$reads)
colnames(no_tad_or_loops) <- c("distance", "reads")
no_tad_or_loops <- ddply(no_tad_or_loops, "distance", summarize, means = mean(reads), sd = sd(reads))

p <- hist(no_tad_or_loops$reads,ylim=range(0,500), breaks=300, plot=TRUE)
p <- hist(no_test_dist_vs_counts_tads$reads, breaks=100, xlim=c(600,1000))
plot(c(2,2,2,2), c(1,2,3,4))
(dist<-fitdistr(as.integer(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)]),densfun = "poisson"))
BIC(dist)
#Poisson fit test
lamb <- trunc(mean(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)]))
hist(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)])
tab1<- table(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)])
test <- rpois(5000, 852)
t <- hist(test)

real <- hist(trunc_1, breaks = 50)
xhist <- c(min(real$breaks), real$breaks)
yhist <- c(min(t$breaks), real$breaks)
xfit<- seq(min(trunc_1), max(trunc_1), length = 40)
yfit<- dpois(xfit,lambda = lamb)
plot(xhist,yhist, type = "s")

freq_exp <- (dpois(0:max(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)]), lambda=lamb))
freq_exp
gf <- goodfit(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)], type = "poisson", method = "MinChisq")
summary(gf)

qqplot(no_tad_or_loops$reads[no_tad_or_loops$distance == (1)])
hist(no_test_dist_vs_counts_tads$reads[no_test_dist_vs_counts_tads$distance == 3])

data_used_vec <- unique(data_used_vec)
graph_pot_tads <- ddply(no_test_dist_vs_counts_tads, "distance", summarize, means = mean(reads), sd = sd(reads))
graph_pot_flares_and_loops <- ddply(flares_and_loops_dvc,"distance", summarize, means = mean(reads), sd = sd(reads))

#plot the signal vs distance curve
plot(dist_vs_avg_counts$distance[dist_vs_avg_counts$distance<50], dist_vs_avg_counts$mean[dist_vs_avg_counts$distance<50],
     xlab = "Distance in 10KB", ylab = "Average Signal", col ="blue", type = "o", main="Signal vs. Distance")
lines(graph_pot_tads$distance[graph_pot_tads$distance < 50], graph_pot_tads$mean[graph_pot_tads$distance < 50], col="red")
lines(graph_pot_flares_and_loops$distance[graph_pot_flares_and_loops$distance < 50], graph_pot_flares_and_loops$mean[graph_pot_flares_and_loops$distance<50], col="green")
lines(no_tad_or_loops$distance[no_tad_or_loops$distance < 50], no_tad_or_loops$means[no_tad_or_loops$distance < 50], col="purple")
legend(20, 1500,c("Global", "Within TADS", "Flare & Loops"), lty = c(1,1), lwd = c(2,2), col=c("blue", "red", "green"))


#### linear regression model ####
tail(dist_vs_avg_counts)
m<- lm( mean ~ distance, data = combined_graph[combined_graph$distance < 50 & combined_graph$model == "Loop&FL",c("distance", "means")] )
summary(m)

m <- nls( mean ~ a/(distance+1), data = dist_vs_avg_counts[dist_vs_avg_counts$distance<20,], start=list(a=2000))
summary(m)

ggplot( dist_vs_avg_counts[dist_vs_avg_counts$distance<50,], aes( x = distance, y= mean ))+
  geom_point()+
  geom_smooth(method = "nls", formula = y ~ a^exp(1)/(x+1),se = F, method.args = list(start = c(a = 3000)))+
  ggtitle("Global Distance vs Mean NLM")+
  theme(plot.title = element_text(hjust = .5))

#combine the graph plots together to make model for all three signal vs distance
no_tad_or_loops$model <- "Background"
graph_pot_flares_and_loops$model <- "Loop&FL"
graph_pot_tads$model <- "TADs"

combined_graph <- rbind(no_tad_or_loops, graph_pot_flares_and_loops, graph_pot_tads)

ggplot( combined_graph[combined_graph$distance < 50,], aes( x = distance, y = means, col=model))+
  geom_line()+
  #geom_smooth(aes(col=model), method = "nls", formula = y ~ a/(x+1), se = F, method.args = list(start = c(a=10000)) )+
  ggtitle("Distance Vs Mean")+
  theme(plot.title = element_text(hjust = .5))


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
  list_of_df <- create_distance_vs_mean_matrices_for_comparison(read_counts, start, end, bot_and_top_vecs)
}


#### Create Location Vectors #####
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


#### Find counts from data ####
create_distance_vs_mean_matrices_for_comparison <- function(reads_mat, start, end, vec_list) {
  d <- end - start
  
  #create matrices
  out_mat <- matrix(ncol = 2)
  in_mat <- matrix(ncol = 2)
  out_mat <- as.data.frame(out_mat)
  in_mat <- as.data.frame(in_mat)
  colnames(out_mat) <- c("distance", "reads")
  colnames(in_mat) <- c("distance", "reads")
  
  condensed_mat <- reads_mat[reads_mat$start >= (start - d) & reads_mat <= (end + d),]
  for ( i in 1:length(condensed_mat) ) {
    #create string that might be found in our condensed matrix
    line_string <- paste(condensed_mat$start[i], ",", condensed_mat$end[i])
    #check to see if it goes in outside matrix
    if ( line_string %in% vec_list[[1]] | line_string %in% vec_list[[2]]) {
      mm_df <- data.frame( (condensed_mat$end[i] - condensed_mat$start[i]),
                           condensed_mat$reads[i] )
      colnames(mm_df) <- colnames(out_mat)
      out_mat <- rbind(out_mat, mm_df)
    }
    #check to see if the info goes inside the matrix
    else if ( line_string %in% vec_list[[3]] | line_string %in% vec_list[[4]] ) {
      mm_df <- data.frame( (condensed_mat$end[i] - condensed_mat$start[i]),
                           condensed_mat$reads[i] )
      colnames(mm_df) <- colnames(in_mat)
      in_mat <- rbind(in_mat, mm_df)
    }
  }
  
  return(list(out_mat, in_mat) )
}
