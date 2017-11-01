#This script analyzes the difference in tad pixel strength as the loops grow in distance

library(MASS)
library(ggplot2)

bin_things <- function(df,bin_total) {
  combined <- c()
  counter <- 0
  for (i in 0:bin_total) {
    combined <- append(combined, i)
    if (length(df$distance[df$distance==i]) == 0) next
    counter <- counter + length(df$distance[df$distance==i])
    if(counter > 50) {
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
  if (length(combined) > 1 && counter > 0) {
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

fit_single_distribution <- function(tad_df, dist_vals) {
  print("Start")
  tad_df$signal <- tad_df$signal+1
  print(dist_vals)
  dist_matrix <- c()
  print("Looping")
  for( i in 0:dist_vals) {
    t <- as.integer(tad_df$signal[tad_df$distance == i])
    #fit the distributions
    nbin_t <- fitdistr(as.integer(t), "negative binomial")
    gam_t <- fitdistr(as.integer(t), "gamma")
    print("t-done")
    #Create data frame from the distributions generated
    
    #For all three at once
    dist_matrix <- append(dist_matrix, c(i,tad_df$model[1],nbin_t$estimate[1], nbin_t$estimate[2], gam_t$estimate[1], gam_t$estimate[2]))
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

create_dist_signal <- function(dataframe, by_val) {
  dist_matrix <- matrix(ncol = 6)
  colnames(dist_matrix) <- c("distance", "model", "theta", "mu", "shape", "rate")
  dist_matrix <- as.data.frame(dist_matrix)
  maxi <- max(dataframe$distance)
  order_model <- c()
  
  for( i in seq(0,maxi,by = by_val) ) {
    dataframe_part <- dataframe[dataframe$Domain_dists > i & dataframe$Domain_dists <= (i+by_val),]
    dataframe_part <- data.frame(dataframe_part$distance, dataframe_part$Signal, paste(as.character(i), "-", as.character(i+by_val), sep = ""))
    order_model <- append(order_model, paste(as.character(i), "-", as.character(i+by_val), sep = ""))
    colnames(dataframe_part) <- c("distance", "signal", "model")
    dataframe_part <- bin_things(dataframe_part, max(dataframe_part$distance))
    dataframe_part <- fit_single_distribution(dataframe_part, max(dataframe_part$distance))
    dataframe_part$model <- as.character(paste(as.character(i), "-", as.character(i+by_val), sep = ""))
    dist_matrix <- rbind(dist_matrix, dataframe_part)
  }
  dist_matrix$model <- factor(dist_matrix$model, order_model)
  
  return(na.omit(dist_matrix))
 }

print_out_data <- function(name, dataframe, o) {
  write.table(dataframe, file = gsub(" ","",paste(o, "/", name, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE, quote = FALSE)
}

real_data <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/count_matrix_chr20_norm.txt", header = TRUE)
real_data$distance <- abs(real_data$Bin2 - real_data$Bin1)/10000
real_data <- na.omit(real_data)
real_data$Domain_dists <- as.character(real_data$Domain_dists)

#have to fix flares
sum(real_data$Flare_1_in > 0)
sum(real_data$Flare_1_out > 0)
sum(real_data$Flare_2_out > 0)
sum(real_data$Flare_2_in > 0)
sum(real_data$Flare_0 > 0)


real_tad <- real_data[real_data$Domains > 0,]
real_flare <- real_data[real_data$Flare_0 > 0 | real_data$Flare_1_in > 0 | real_data$Flare_2_in > 0 | real_data$Flare_1_out > 0 | real_data$Flare_2_out > 0, ]
real_bg <- real_data[real_data$Flare_0 == 0 & real_data$Flare_1_in == 0 & real_data$Flare_2_in == 0 &
                       real_data$Flare_1_out == 0 & real_data$Flare_2_out == 0 & real_data$Domains == 0 & real_data$Loops == 0, ]

tad1 <- real_data[real_data$Domains == 1,]
tad1$Domain_dists <- as.integer(gsub(",","",tad1$Domain_dists))

small <- tad1[tad1$Domain_dists <= 20,]
large <- tad1[tad1$Domain_dists > 20 & tad1$Domain_dists <= 40,]
reg <- tad1[tad1$Domain_dists>40 & tad1$Domain_dists <= 60,]

small <- data.frame(small$distance, small$Signal, "TADs")
large <- data.frame(large$distance, large$Signal, "Loop&FL")
reg <- data.frame(reg$distance, reg$Signal, "Background")

#get correct names of real data 
colnames(small) <- c("distance", "signal", "model")
colnames(large) <- c("distance", "signal", "model")
colnames(reg) <- c("distance", "signal", "model")

small$model <- as.character(small$model)
large$model <- as.character(large$model)
reg$model <- as.character(reg$model)


small <- bin_things(small,20)
large <- bin_things(large,20)
reg <- bin_things(reg,20)
fit <- fit_single_distribution(small,20)
fit2 <- fit_single_distribution(large, 40)
fit3 <- fit_single_distribution(reg,60)
fit$model[fit$model == "TADs"] <- "20"
fit2$model[fit2$model == "TADs"] <- "20-40"
fit3$model[fit3$model == "TADs"] <- "40-60"
ggplot(fit, aes(x=distance, y=mu, col=model))+
  geom_line()+
  geom_line(data=fit2, aes(x=distance, y=mu, col=model))+
  geom_line(data=fit3, aes(x=distance, y=mu, col=model))


#### This code is to create a distance to signal for diff loop values ####
tad1 <- real_data[real_data$Domains == 1,]
tad1$Domain_dists <- as.integer(gsub(",","",tad1$Domain_dists))
new_fit <- create_dist_signal(tad1,20)

ggplot(new_fit, aes(x=distance, y=mu, col=model))+
  geom_line()+
  ggtitle("TAD Signal Decay is Dependant on Loop Size")+
  theme( plot.title = element_text(hjust = .5) )+
  labs( x="Linear Genomic Distance (Mb)", y= "Normalized Mean Counts", col = "Loop Distance (Mb)")+
  scale_y_continuous(limits = c(0,300))


outer <- "/Users/ncolaian/Documents/phanstiel_lab/data/"
print_out_data("tad_sig_loop_dist", new_fit, outer)


#This code is for looking at the flares
real_flare <- real_data[real_data$Flare_0 > 0 | real_data$Flare_1_in > 0 | real_data$Flare_2_in > 0 | real_data$Flare_1_out > 0 | real_data$Flare_2_out > 0, ]
real_flare <- real_flare[,c("Signal", "distance", "Flare_0", "Flare_1_in", "Flare_1_out", "Flare_2_in", "Flare_2_out", "Domain_dists")]
#just look at 0 flares
zero_flare <- real_data[real_data$Domains ==1 & real_data$Flare_0 == 1 & real_data$Flare_1_in == 0 & real_data$Flare_2_in == 0 & real_data$Flare_1_out == 0 & real_data$Flare_2_out == 0, ]

real_flare <- melt(real_flare, id.vars = c("Signal", "distance","Domain_dists"))
real_flare$variable <- as.character(real_flare$variable)
real_flare$value[real_flare$value == 0] <- NA
real_flare <- na.omit(real_flare)

zero_flare$Domain_dists <- as.integer(gsub(",","",zero_flare$Domain_dists))
rf_0 <- create_dist_signal(zero_flare, 40)


ggplot(real_flare, aes(distance, log(Signal), col = variable)) +
  geom_line()

ggplot(rf_0, aes(x=distance, y=mu, col=model))+
  geom_line()+
  ggtitle("Flares 0 Away Signal vs Distance")+
  theme( plot.title = element_text(hjust = .5) )+
  labs( x="Pixel Distance (Mb)", y= "Signal", col = "Loop Distance")
  scale_y_continuous(limits = c(0,300))
