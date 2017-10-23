#This script analyzes the difference in tad pixel strength as the loops grow in distance

library(MASS)
library(ggplot2)
library(plyr)
library(reshape2)

bin_things <- function(df,bin_total) {
  combined <- c()
  counter <- 0
  for (i in seq(0:bin_total)) {
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
            if (nrow(mm) == 0) next
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
        if (nrow(mm) == 0) next
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
    print(i)
    if ( NA %in% t || length(t) <= 1) next
    #fit the distributions
    nbin_t <- fitdistr(as.integer(t), "negative binomial")
    gam_t <- fitdistr(as.integer(t), "gamma")
    print("t-done")
    #Create data frame from the distributions generated
    
    #For all three at once
    dist_matrix <- append(dist_matrix, c(i,tad_df$model[1],nbin_t$estimate[1], nbin_t$estimate[2], gam_t$estimate[1], gam_t$estimate[2]))
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
  dataframe <- arrange(dataframe, loop_sig)
  uq_ls <- unique(dataframe$loop_sig)
  
  last <- 1
  for( i in seq(as.integer(length(uq_ls)/by_val),length(uq_ls),by = as.integer(length(uq_ls)/by_val)) ) {
    sigs <- uq_ls[last:i]
    dataframe_part <- dataframe[dataframe$loop_sig %in% sigs,]
    distance <- ((dataframe_part$Bin2 - dataframe_part$Bin1)/10000)
    if (nrow(dataframe_part) <= 1) next
    maxi <- as.integer(max(dataframe_part$loop_sig))
    mini <- as.integer(min(dataframe_part$loop_sig))
    dataframe_part <- data.frame(distance, dataframe_part$Signal, as.character(maxi))
    colnames(dataframe_part) <- c("distance", "signal", "model")
    dataframe_part <- bin_things(na.omit(dataframe_part), max(dataframe_part$distance))
    dataframe_part <- fit_single_distribution(dataframe_part, max(dataframe_part$distance))
    if ( nrow(dataframe_part) == 0) next
    dataframe_part$model <- paste(as.character(mini), "-", as.character(maxi), sep = "")
    dist_matrix <- rbind(dist_matrix, dataframe_part)
    last <- i+1
  }
  
  return(na.omit(dist_matrix))
}

print_out_data <- function(name, dataframe, o) {
  write.table(dataframe, file = gsub(" ","",paste(o, "/", name, ".txt",sep = ""),fixed = TRUE), sep = "\t", row.names = FALSE, quote = FALSE)
}

get_areas_in_tads_flares_loops <- function(need_loop, real, bl,flares=TRUE, tads=TRUE) {
  loop_vec <- c()
  flare_vec <- c()
  tad_vec <- c()
  tad_loop <- c()
  flare_loop <- c()
  
  for(i in seq(1:length(need_loop$Bin1))) {
    starting <- need_loop$Bin1[i]
    ending <- need_loop$Bin2[i]
    dist <- (ending-starting)/bl
    sig <- real$Signal[real$Bin1 == starting & real$Bin2 == ending]
    
    #always get loops
    loop_vec <- append(loop_vec, paste(starting, ":", ending, sep = ""))
    
    #get both flares and tads
    if(flares == TRUE && tads == TRUE) {
      #f_val and t_val keep track of how many tad pixels and flare pixels are added
      f_val <- length(flare_vec)
      t_val <- length(tad_vec)
      flare_vec <- append(flare_vec, get_flares(starting, ending, bl))
      tad_vec <- append(tad_vec, get_tads(starting,ending, bl))
      f_val <- length(flare_vec) - f_val
      t_val <- length(tad_vec) - t_val
      tad_loop <- append(tad_loop, rep(sig,t_val))
      flare_loop <- append(flare_loop, rep(sig,f_val))
    }
    #just get tads
    else if( tads == TRUE && flares == FALSE ) {
      t_val <- length(tad_vec)
      tad_vec <- append(tad_vec, get_tads(starting,ending, bl, fl = FALSE))
      t_val <- length(tad_vec) - t_val
      tad_loop <- append(tad_loop, rep(sig,t_val))
    }
    #just get flares
    else{
      f_val <- length(flare_vec)
      flare_vec <- append(flare_vec, get_flares(starting, ending, bl))
      f_val <- length(flare_vec) - f_val
      flare_loop <- append(flare_loop, rep(sig,f_val))
    }
  }
  
  return(list(loop_vec, flare_vec, tad_vec, flare_loop, tad_loop))
}

create_combined_tag <- function(x,spot) {
  tag <- paste(spot,":",x, sep="")
  return(tag)
}

get_flares <- function(starting, ending, bl) {
  vec<- c()
  adjust <- bl
  flare_end <- ending+adjust
  flare_start <- starting - adjust
  bot <- seq(from = flare_start, to = abs(starting+adjust), by = bl)
  top <- seq(to = flare_end, from = abs(ending-adjust), by = bl)
  
  #get left flare and top right combined
  for(spot in bot) {
    new <- sapply(seq(spot, flare_end, by = bl), function(x) {x <- paste(spot,":",x, sep=""); return(x)})
    vec <- append(vec, new)
  }
  
  #get top flare
  for(spot in top) {
    new <- sapply(seq(to = spot, from = starting, by = bl), function(x) {x <- paste(x,":",spot, sep=""); return(x)})
    vec <- append(vec, new)
  }
  
  return(vec)
}

get_tads <- function(starting, ending, bl, fl=TRUE) {
  vec <- c()
  new_s <- starting
  new_e <- ending
  if (fl == TRUE) {
    adjust <- bl
    new_s <- starting+adjust
    new_e <- ending-adjust
  }
  
  for(spot in seq(from = new_s, to = new_e, by=bl)) {
    new <- sapply(seq(spot, new_e, by = bl), function(x) {x <- paste(spot,":",x, sep=""); return(x)})
    vec <- append(vec, new)
  }
  
  return(vec)
}
#file on personal laptop
#real_data <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/count_matrix_chr20_old.txt", header = TRUE)
#file on shared computer
real_data <- read.delim("/Users/phanstiel4/Documents/code_rep/data/count_matrix_chr20.txt", header = TRUE)

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
loops_to_sim <- real_data[real_data$Loops == 1,]
non_bg <- get_areas_in_tads_flares_loops(loops_to_sim,real_data, 10000, flares = FALSE, tads = TRUE)

tad1 <- na.omit(real_data[real_data$Domains == 1 & create_combined_tag(real_data$Bin2, real_data$Bin1) %in% non_bg[[3]],])
tad1$loop_sig <- non_bg[[5]][which(non_bg[[3]] %in% paste(tad1$Bin1,":", tad1$Bin2, sep = ""))]
by <- as.integer((max(tad1$loop_sig) - min(tad1$loop_sig))/10)
tad1$loop_sig <- as.integer(tad1$loop_sig)
new_fit <- create_dist_signal(tad1,5)

ggplot(new_fit, aes(x=distance, y=log(mu), col=model))+
  geom_line()+
  ggtitle("Tad Signal vs Distance")+
  theme( plot.title = element_text(hjust = .5) )+
  labs( x="Distance (Mb)", y= "Log(Loop Signal)", col = "Loop Distance")


outer <- "/Users/ncolaian/Documents/phanstiel_lab/data/"
print_out_data("tad_sig_loop_dist", new_fit, outer)





