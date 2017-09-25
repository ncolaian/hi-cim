#! /usr/bin/env Rscript

library(MASS)
library(lmom)

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

#### Subroutines ####
fit_distributions <- function(tad_df, loop_df, back_df) {
  tad_df <- na.omit(tad_df)
  loop_df <- na.omit(loop_df)
  back_df <- na.omit(back_df)
  
  dist_vals <- seq(1,30)
  dist_matrix <- c()
  
  for( i in dist_vals) {
    t <- as.integer(tad_df$reads[tad_df$distance == i])
    lf <- as.integer(loop_df$reads[loop_df$distance == i])
    b <- as.integer(back_df$reads[back_df$distance == i])
    print(t)
     #fit the distributions
    nbin_t <- fitdistr(as.integer(t), "negative binomial")
    gam_t <- fitdistr(as.integer(t), "gamma")
    
    nbin_lf <- fitdistr(as.integer(lf), "negative binomial")
    gam_lf <- fitdistr(as.integer(t), "gamma")
    
    b[is.na(b)] <- 0
    b <- b+1
    nbin_b <- fitdistr(as.integer(b), "negative binomial")
    gam_b <- fitdistr(as.integer(b), "gamma")
    
    #Create data frame from the distributions generated
    
    #For all three at once
    dist_matrix <- append(dist_matrix, c(i,"TADs",nbin_t$estimate[1], nbin_t$estimate[2], gam_t$estimate[1], gam_t$estimate[2],
                                         i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2],
                                         i, "Background", nbin_b$estimate[1], nbin_b$estimate[2], gam_b$estimate[1], gam_b$estimate[2]) )
    
    #for a single one at a time
    #dist_matrix <- append(dist_matrix, c(i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2]))
  }
  
  #create a matrix from the vector
  dist_matrix <- matrix(dist_matrix, ncol = 6, byrow = TRUE)
  colnames(dist_matrix) <- c("distance", "model", "theta", "mu", "scale", "rate")
  
  #correct variable problems
  dist_matrix <- as.data.frame(dist_matrix)
  dist_matrix$distance <- as.integer(dist_matrix$distance)
  dist_matrix$model <- as.character(dist_matrix$model)
  dist_matrix$theta <- as.numeric(dist_matrix$theta)
  dist_matrix$mu <- as.numeric(dist_matrix$mu)
  dist_matrix$scale <- as.numeric(dist_matrix$scale)
  dist_matrix$rate <- as.numeric(dist_matrix$rate)
  
  return(as.data.frame(dist_matrix))
}

print_out_data <- function(name, dataframe, o) {
  write.table(dataframe, file = gsub(" ","",paste(o, "/", name, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE, quote = FALSE)
}


#### MAIN ####

#the combined data points
data4tad <- read.delim(args[1], quote = "")
data4loop <- read.delim(args[2], quote = "")
data4back <- read.delim(args[3], quote = "")
#testing
#data4distr <- read.delim("/Users/phanstiel4/Documents/sim_graphs/distribution/comb2_d.txt", stringsAsFactors = FALSE)
#get output file path
outer <- args[4]
colnames(data4tad) <- c("distance", "reads", "used", "model")
colnames(data4loop) <- c("distance", "reads", "used", "model")
colnames(data4back) <- c("distance", "reads", "used", "model")

data4tad$distance <- as.integer(data4tad$distance)
data4loop$distance <- as.integer(data4loop$distance)
data4back$distance <- as.integer(data4back$distance)

#make dataframe
fit <- fit_distributions(data4tad, data4loop, data4back)

#print out table
print_out_data("distr_50", fit, outer)
#test
#outer <- "/Users/phanstiel4/Documents/sim_graphs/distribution"
#print_out_data("distr_50", fit, outer)
