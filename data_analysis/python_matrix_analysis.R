#This script will create some distribution files and also graphs to look at the data collected
#using the python matrix builder

library(MASS)
library(ggplot2)
library(plyr)

#### Subroutines ####
fit_distributions <- function(tad_df, loop_df, back_df) {
  print("Start")
  tad_df$signal <- tad_df$signal+1
  loop_df$signal <- loop_df$signal+1
  back_df$signal <- back_df$signal+1
  
  dist_vals <- 192
  dist_matrix <- c()
  print("Looping")
  for( i in 0:dist_vals) {
    t <- as.integer(tad_df$signal[tad_df$distance == i])
    lf <- as.integer(loop_df$signal[loop_df$distance == i])
    b <- as.integer(back_df$signal[back_df$distance == i])
    
    
    
    #fit the distributions
    nbin_t <- fitdistr(as.integer(t), "negative binomial")
    gam_t <- fitdistr(as.integer(t), "gamma")
    print("t-done")
    nbin_lf <- fitdistr(as.integer(lf), "negative binomial")
    gam_lf <- fitdistr(as.integer(lf), "gamma")
    print("lf_done")
    
    b <- b+1
    nbin_b <- fitdistr(as.integer(b), "negative binomial")
    print("b-bin")
    gam_b <- fitdistr(as.integer(b), "gamma")
    print("b-gam")
    #Create data frame from the distributions generated
    
    #For all three at once
    dist_matrix <- append(dist_matrix, c(i,"TADs",nbin_t$estimate[1], nbin_t$estimate[2], gam_t$estimate[1], gam_t$estimate[2],
                                         i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2],
                                         i, "Background", nbin_b$estimate[1], nbin_b$estimate[2], gam_b$estimate[1], gam_b$estimate[2]) )
    print(i)
    #for a single one at a time
    #dist_matrix <- append(dist_matrix, c(i, "Loop&FL", nbin_lf$estimate[1], nbin_lf$estimate[2], gam_lf$estimate[1], gam_lf$estimate[2]))
  }
  
  #create a matrix from the vector
  dist_matrix <- matrix(dist_matrix, ncol = 6, byrow = TRUE)
  colnames(dist_matrix) <- c("distance", "model", "theta", "mu", "shape", "rate")
  
  #correct variable problems
  dist_matrix <- as.data.frame(dist_matrix)
  dist_matrix$distance <- as.integer(as.character(dist_matrix$distance))
  dist_matrix$model <- as.character(dist_matrix$model)
  dist_matrix$theta <- as.numeric(as.character(dist_matrix$theta))
  dist_matrix$mu <- as.numeric(as.character(dist_matrix$mu))
  dist_matrix$scale <- as.numeric(as.character(dist_matrix$scale))
  dist_matrix$rate <- as.numeric(as.character(dist_matrix$rate))
  
  return(as.data.frame(dist_matrix))
}

print_out_data <- function(name, dataframe, o) {
  write.table(dataframe, file = gsub(" ","",paste(o, "/", name, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE, quote = FALSE)
}

bin_things <- function(df,bin_total) {
  combined <- c()
  counter <- 0
  for (i in 0:bin_total) {
    combined <- append(combined, i)
    if (length(df$distance[df$distance==i]) == 0) next
    counter <- counter + length(df$distance[df$distance==i])
    if(counter > 450) {
      if (length(combined) > 1) {
        add <- matrix(ncol = ncol(df))
        colnames(add) <- colnames(df)
        for (j in combined) {
          for (l in combined) {
            if (l == j) next
            mm <- df[df$distance == j, ]
            mm$distance <- l
            add <- rbind(add,mm)
          }
        }
        df <- rbind(df, na.omit(add))
      }
      combined <- c()
      counter <- 0
    }
  }
  #combine left over bins
  if (length(combined) > 1) {
    add <- matrix(ncol = ncol(df))
    colnames(add) <- colnames(df)
    for (j in combined) {
      for (l in combined) {
        if (l == j) next
        mm <- df[df$distance == j, ]
        mm$distance <- l
        add <- rbind(add,mm)
      }
    }
    df <- rbind(df, na.omit(add))
  }
  
  return(df)
}



#WORK WITH NEW DATAFRAME
real_data <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/count_matrix_chr20.txt", header = TRUE)
real_data$distance <- abs(real_data$Bin2 - real_data$Bin1)/10000
real_data <- na.omit(real_data)

real_tad <- real_data[real_data$Domains > 0,]
real_flare <- real_data[real_data$Flare_0 > 0 | real_data$Flare_1 > 0 | real_data$Flare_2 > 0, ]
real_bg <- real_data[real_data$Flare_0 == 0 & real_data$Flare_1 == 0 & real_data$Flare_2 == 0 &
                       real_data$Domains == 0 & real_data$Loops == 0, ]

real_tad <- data.frame(real_tad$distance, real_tad$Signal, "TADs")
real_flare <- data.frame(real_flare$distance, real_flare$Signal, "Loop&FL")
real_bg <- data.frame(real_bg$distance, real_bg$Signal, "Background")

#get correct names of real data 
colnames(real_tad) <- c("distance", "signal", "model")
colnames(real_flare) <- c("distance", "signal", "model")
colnames(real_bg) <- c("distance", "signal", "model")

real_tad$model <- as.character(real_tad$model)
real_flare$model <- as.character(real_flare$model)
real_bg$model <- as.character(real_bg$model)

#real_tad <- ddply(real_tad, .(distance,model), summarize, pixels = NROW(signal) ,mean_sig = mean(signal))
#real_flare <- ddply(real_flare, .(distance,model), summarize, pixels = NROW(signal) ,mean_sig = mean(signal))
#real_bg <- ddply(real_bg, .(distance,model), summarize, pixels = NROW(signal) ,mean_sig = mean(signal))

#bin real data
real_tad <- bin_things(real_tad, 192)
real_flare <- bin_things(real_flare,192)
real_bg <- bin_things(real_bg,192)

fit <- fit_distributions(real_tad, real_flare, real_bg)

max(real_bg$distance)
ggplot(fit, aes(x=distance, y=scale/rate^2, col=model))+
  geom_line()
hist(real_bg$signal[real_bg$distance==150])

outer <- "/Users/ncolaian/Documents/phanstiel_lab/data/"
print_out_data("distr_192", fit, outer)
