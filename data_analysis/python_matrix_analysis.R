#This script will create some distribution files and also graphs to look at the data collected
#using the python matrix builder

library(MASS)
library(ggplot2)
library(plyr)

#### Subroutines ####
fit_distributions <- function(tad_df, loop_df, back_df, dist_vals) {
  print("Start")
  tad_df$signal <- tad_df$signal+1
  loop_df$signal <- loop_df$signal+1
  back_df$signal <- back_df$signal+1
  
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
  dist_matrix$shape <- as.numeric(as.character(dist_matrix$shape))
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

get_log_fold_changes <- function(means_df) {
  types <- unique(means_df$model)
  return_df_distance <- c()
  return_df_model <- c()
  return_df_fold_change <- c()
  
  #compares each type of moel to each other at each distance
  for(i in seq(1,length(types)) ) {
    for ( j in seq(1,length(types)) ) {
      if ( i == j ) next # Handles comparing it to itself
      #find max distance that both models have values for
      maxi <- min(c(max(means_df$distance[means_df$model == types[i]]), max(means_df$distance[means_df$model == types[j]])) )
      
      #go through and get fold change for each distance
      for( d in seq(0,maxi) ) {
        fc <- ( (means_df$mu[means_df$model == types[i] & means_df$distance == d])/(means_df$mu[means_df$model == types[j] & means_df$distance == d]) )
        return_df_distance <- append(return_df_distance, d)
        return_df_model <- append(return_df_model, paste(types[i], "/", types[j], sep = ""))
        return_df_fold_change <- append(return_df_fold_change, as.numeric(fc))
      }
    }
  }
  return_df <- data.frame(return_df_distance, return_df_model, return_df_fold_change)
  colnames(return_df) <- c("distance", "model", "fold_change")
  return_df$distance <- as.integer(return_df$distance)
  return_df$model <- as.character(return_df$model)
  return_df$fold_change <- as.numeric(return_df$fold_change)
  return( na.omit(return_df) )
}



#WORK WITH NEW DATAFRAME
real_data <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/count_matrix_chr20_norm.txt", header = TRUE)
real_data$distance <- abs(real_data$Bin2 - real_data$Bin1)/10000
real_data <- na.omit(real_data)

real_tad <- real_data[real_data$Domains > 0,]
real_flare <- real_data[(real_data$Flare_0 > 0 | real_data$Flare_1_in > 0 | real_data$Flare_2_in > 0 | real_data$Flare_1_out > 0 | real_data$Flare_2_out > 0) & real_data$Loops == 0, ]
real_bg <- real_data[real_data$Flare_0 == 0 & real_data$Flare_1_in == 0 & real_data$Flare_2_in == 0 &
                       real_data$Flare_1_out == 0 & real_data$Flare_2_out == 0 & real_data$Domains == 0 & real_data$Loops == 0, ]
real_loop <- real_data[real_data$Loops > 0,]

real_tad <- data.frame(real_tad$distance, real_tad$Signal, "TADs")
real_flare <- data.frame(real_flare$distance, real_flare$Signal, "Loop&FL")
real_bg <- data.frame(real_bg$distance, real_bg$Signal, "Background")
real_loop <- data.frame(real_loop$distance,"Loop", real_loop$Signal)

#get correct names of real data 
colnames(real_tad) <- c("distance", "signal", "model")
colnames(real_flare) <- c("distance", "signal", "model")
colnames(real_bg) <- c("distance", "signal", "model")
colnames(real_loop) <- c("distance", "model", "mu")

real_tad$model <- as.character(real_tad$model)
real_flare$model <- as.character(real_flare$model)
real_bg$model <- as.character(real_bg$model)
real_loop$model <- as.character(real_loop$model)

#real_tad <- ddply(real_tad, .(distance,model), summarize, pixels = NROW(signal) ,mean_sig = mean(signal))
#real_flare <- ddply(real_flare, .(distance,model), summarize, pixels = NROW(signal) ,mean_sig = mean(signal))
#real_bg <- ddply(real_bg, .(distance,model), summarize, pixels = NROW(signal) ,mean_sig = mean(signal))

#bin real data
real_tad <- bin_things(real_tad, 191)
real_flare <- bin_things(real_flare,191)
real_bg <- bin_things(real_bg,191)

fit <- fit_distributions(real_tad, real_flare, real_bg, 191)

#Plot the means
fit <- data.frame(fit$distance, fit$model, fit$mu)
colnames(fit) <- c("distance", "model", "mu")
fit$model <- as.character(fit$model)
fit$model[fit$model == "Loop&FL"] <- "Flares"

#combine the loops
real_loop <- ddply(real_loop, .(distance,model), summarise, mu = mean(mu))
fit <- rbind(fit,real_loop)

ggplot(fit, aes(x=distance, y=log(mu), col=model))+
  geom_line()+
  ggtitle("Normalized Counts Decrease at Different Rates")+
  labs( x="Linear Genomic Distance (Mb)", y= "Log(Mean Counts)", color="Feature")+
  theme( plot.title = element_text(hjust = .5),
       legend.title = element_text(hjust = .5))
  
  
hist(real_bg$signal[real_bg$distance==150])

outer <- "/Users/ncolaian/Documents/phanstiel_lab/data/"
print_out_data("distr_192", fit, outer)

#### TAD ANALYSIS ####
#This is analysis of the different tads
tad1 <- real_data[real_data$Domains == 1,]
tad2 <- real_data[real_data$Domains == 2,]
tad_great <- real_data[real_data$Domains > 2,]


tad1 <- data.frame(tad1$distance, tad1$Signal, "TADs")
tad2 <- data.frame(tad2$distance, tad2$Signal, "Loop&FL")
tad_great <- data.frame(tad_great$distance, tad_great$Signal, "Background")

#get correct names of real data 
colnames(tad1) <- c("distance", "signal", "model")
colnames(tad2) <- c("distance", "signal", "model")
colnames(tad_great) <- c("distance", "signal", "model")

tad1$model <- as.character(tad1$model)
tad2$model <- as.character(tad2$model)
tad_great$model <- as.character(tad_great$model)

#analysis
tad1_b <- bin_things(tad1, 146)
tad2_b <- bin_things(tad2,146)
tad_great_b <- bin_things(tad_great,146)
fit <- fit_distributions(tad1_b, tad2_b, tad_great_b, 146)
fit$model[fit$model == "Loop&FL"] <- "2 TADs"
fit$model[fit$model == "Background"] <- "3+ TADs"

#fold change analysis
fold_changes <- get_log_fold_changes(fit)

#subset the fold changes I want
wanted_fold_changes <- fold_changes[fold_changes$model == "2 TADs/TADs" | fold_changes$model == "3+ TADs/2 TADs",]

#plot original fit models to investigate mean pix values over distance with differing tad amounts
ggplot(fit, aes(x=distance, y=log(mu), col=model))+
  geom_line()+
  ggtitle("Normalized Mean Counts Per Pixel vs Distance")+
  theme( plot.title = element_text(hjust = .5) )+
  labs( x="Distance (Mb)", y= "Log(Mean)", color="Number of TADs")

#Plot the fold changes over distance
ggplot(wanted_fold_changes, aes(x=distance, y=fold_change, col=model))+
  geom_line()
  

#Checking the point to cut off small vs large
fit1 <- fit[fit$model == "TADs",]
fit2 <- fit[fit$model == "TAD2",]
for ( i in unique(fit$distance) ) {
  variance <- (fit1$shape[fit1$distance == i]/(fit1$rate[fit1$distance == i]^2))
  sd <- sqrt(variance)
  sd <- log(sd)
  print(c(log(fit1$mu[fit1$distance == i]),log(fit2$mu[fit2$distance == i])))
  if ( log(fit1$mu[fit1$distance == i])+sd < log(fit2$mu[fit2$distance == i])) {
    break
  }
}

outer <- "/Users/ncolaian/Documents/phanstiel_lab/data/"
print_out_data("tad_distr_146", fit, outer)
