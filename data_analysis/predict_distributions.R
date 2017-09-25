#! /usr/bin/env Rscript

library(MASS)
library(lmom)

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

#### Subroutines ####
fit_distributions <- function(comb_df) {
  comb_df <- na.omit(comb_df)
  dist_vals <- seq(1,30)
  dist_matrix <- c()
  
  for( i in dist_vals) {
    #t <- as.integer(comb_df$reads[comb_df$model == "TADs" & comb_df$distance == i])
    #lf <- as.integer(comb_df$reads[comb_df$model == "Loop&FL" & comb_df$distance == i])
    #b <- as.integer(comb_df$reads[comb_df$model == "Background" & comb_df$distance == i])
    # 
    # #fit the distributions
    # nbin_t <- fitdistr(as.integer(t), "negative binomial")
    # gam_t <- fitdistr(as.integer(t), "gamma")
    
    nbin_lf <- fitdistr(as.integer(comb_df$reads[comb_df$distance == i]), "negative binomial")
    gam_lf <- fitdistr(as.integer(comb_df$reads[comb_df$distance == i]), "gamma")
    
    # b[is.na(b)] <- 0
    # b <- b+1
    # nbin_b <- fitdistr(as.integer(b), "negative binomial")
    # gam_b <- fitdistr(as.integer(b), "gamma")
    
    #Create data frame from the distributions generated
    
    #For all three at once
    #dist_matrix <- append(dist_matrix, c(i,"TADs",nbin_t$estimate[1], nbin_t$estimate[2], gam_t$estimate[1], gam_t$estimate[2],
                                         #i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2],
                                        # i, "Background", nbin_b$estimate[1], nbin_b$estimate[2], gam_b$estimate[1], gam_b$estimate[2]) )
    
    #for a single one at a time
    dist_matrix <- append(dist_matrix, c(i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2]))
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
data4distr <- read.delim(args[1], quote = "")
#testing
#data4distr <- read.delim("/Users/phanstiel4/Documents/sim_graphs/distribution/comb2_d.txt", stringsAsFactors = FALSE)
#get output file path
outer <- args[2]
colnames(data4distr) <- c("distance", "reads", "used", "model")
data4distr$distance <- as.integer(data4distr$distance)

#make dataframe
fit <- fit_distributions(data4distr)

#print out table
print_out_data("distr_50", fit, outer)
#test
#outer <- "/Users/phanstiel4/Documents/sim_graphs/distribution"
#print_out_data("distr_50", fit, outer)
