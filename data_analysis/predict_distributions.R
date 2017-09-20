#! /usr/bin/env Rscript

library(MASS)

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

#### Subroutines ####
fit_distributions <- function(comb_df) {
  dist_vals <- seq(1,50)
  dist_matrix <- c()
  
  for( i in dist_vals) {
    t <- comb_df$reads[comb_df$model == "TADs" & comb_df$distance == i]
    lf <- comb_df$reads[comb_df$model == "Loop&FL" & comb_df$distance == i]
    b <- comb_df$reads[comb_df$model == "Background" & comb_df$distance == i]
    colnames(t) <- c("reads")
    colnames(lf) <- c("reads")
    colnames(b) <- c("reads")
    
    #fit the distributions
    nbin_t <- fitdistr(na.omit(as.integer(t)), "negative binomial")
    gam_t <- fitdistr(na.omit(as.integer(t)), "gamma")
    
    nbin_lf <- fitdistr(na.omit(as.integer(lf)), "negative binomial")
    gam_lf <- fitdistr(na.omit(as.integer(lf)), "gamma")
    
    nbin_b <- fitdistr(na.omit(as.integer(b)), "negative binomial")
    gam_b <- fitdistr(na.omit(as.integer(b)), "gamma")
    
    #Create data frame from the distributions generated
    nb_t_df <- rnegbin(length(t$reads), mu = nbin_t$estimate[2], theta = nbin_t$estimate[1])
    gam_t_df <- rgamma(length(t$reads), gam_t$estimate[1], rate = gam_t$estimate[2])
    
    dist_matrix <- append(dist_matrix, c(i,"TADs",nbin_t$estimate[1], nbin_t$estimate[2], gam_t$estimate[1], gam_t$estimate[2],
                                         i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2],
                                         i, "Background", nbin_b$estimate[1], nbin_b$estimate[2], gam_b$estimate[1], gam_b$estimate[2]) )
  }
  
  #create a matrix from the vector
  dist_matrix <- matrix(dist_matrix, ncol = 6, byrow = TRUE)
  colnames(dist_matrix) <- c("distance", "model", "theta", "mu", "scale", "rate")
  return(as.data.frame(dist_matrix))
}

print_out_data <- function(name, dataframe, out) {
  write.table(dataframe, file = gsub(" ","",paste(out, "/", name, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE)
}


#### MAIN ####

#the combined data points
data4distr <- read.delim(args[1])
#get output file path
out <- args[2]

#make dataframe
fit <- fit_distributions(na.omit(data4distr))

#print out table
print_out_data("distr_50", fit, out)