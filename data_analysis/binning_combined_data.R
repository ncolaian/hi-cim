# This script will attempt to spline fit binned data

library(ggplot2)
library(plyr)

comb<- read.delim("/Users/phanstiel4/Documents/sim_graphs/second_pass/full_sig_distance_all_things.txt", header = TRUE)

#You need about 450 reads to get a good bin

find_middle_of_two_bins <- function(N1, N2, end_pix) {
  total <- N1 + N2
  point <- end_pix - (N1/total)
  return(point)
}

find_middle_point_of_mult_bins <- function(bin_counts_df) {
  find_middle_mat <- matrix(ncol = 4)
  colnames(find_middle_mat) <- c("dist", "btot", "etot", "diff")
  find_middle_mat <- as.data.frame(find_middle_mat)
  
  #Goes through and creates dataframe that is used to find two closest means as defined as having simalar pixel counts on each side
  for( i in seq(1:nrow(bin_counts_df) ) ) {
    before_tot <- sum(as.numeric(bin_counts_df[1:i,5]), na.rm = TRUE)
    end_tot <- sum(as.numeric(bin_counts_df[i:nrow(bin_counts_df),5]), na.rm = TRUE)
    mm <- matrix(c(bin_counts_df[i,1], before_tot, end_tot, abs(before_tot-end_tot)), ncol = 4)
    colnames(mm) <- colnames(find_middle_mat)
    find_middle_mat <- rbind(find_middle_mat, mm)
  }
  
  
  find_middle_mat <- na.omit(find_middle_mat)
  minval <- min(find_middle_mat$diff)
  middle_p <- find_middle_mat[find_middle_mat$diff == minval,]
  middle_p <- as.numeric(middle_p)
  distance <- find_middle_of_two_bins(as.numeric(middle_p[2]), as.numeric(middle_p[3]), (as.numeric(middle_p[1])+1))
  mean_of_bins <- sum(as.numeric(bin_counts_df[,2]))/sum(as.numeric(bin_counts_df[,5]))
  vec_values <- c(distance, mean_of_bins)
  
  return(vec_values)
}

combine_bins_if_necessary <- function(min_read_val, df) {
  new_df <- matrix(ncol = 2)
  colnames(new_df) <- c("distance", "means")
  new_df <- as.data.frame(na.omit(new_df))
  
  middle_df <- matrix(ncol = 7)
  middle_df = as.data.frame( middle_df, col.names = colnames(df) )
  middle_df <- na.omit(middle_df)
  
  #go through 100 bins
  for( i in seq(0,100) ) {
    #make sure that distance 
    if ( !(i %in% df$distance) ) {
      mm <- matrix(c(i,0,0,0,0,as.character(df$model[1]),0), ncol = 7)
      colnames(mm) <- colnames(df)
      middle_df <- rbind(middle_df, mm)
      next
    }
    
    middle_df <- rbind(middle_df, df[df$distance == i,])
    #this means there are enough pixel
    if(sum(as.numeric(middle_df[,5]),na.rm = TRUE) >= min_read_val) {
      
      if(nrow(middle_df) == 1) {
        mm <- matrix(middle_df[1,c(1,3)], ncol = 2)
        colnames(mm) <- c("distance", "means")
        new_df <- rbind(new_df, as.data.frame(mm))
      }
      else if( nrow(middle_df) == 2 ) {
        distance <- find_middle_of_two_bins(middle_df[1,5], middle_df[2,5], i)
        m <- sum(middle_df[,2])/sum(middle_df[,5])
        mm <- matrix(c(distance,m), ncol = 2)
        colnames(mm) <- c("distance", "means")
        new_df <- rbind(new_df, as.data.frame(mm))
      }
      else{
        vec <- find_middle_point_of_mult_bins(middle_df)
        mm <- matrix(vec, ncol = 2)
        colnames(mm) <- c("distance", "means")
        new_df <- rbind(new_df, as.data.frame(mm))
      }
      
      middle_df <- middle_df[0,] #reset dataframe holding previous records
    } #enough pixels if
    
    #handle last thing
    if( i == 100 && nrow(middle_df) > 0) {
      if(nrow(middle_df) == 1) {
        mm <- matrix(middle_df[1,c(1,3)], ncol = 2)
        colnames(mm) <- c("distance", "means")
        new_df <- rbind(new_df, as.data.frame(mm))
      }
      else{
        vec <- find_middle_point_of_mult_bins(middle_df)
        mm <- matrix(vec, ncol = 2)
        colnames(mm) <- c("distance", "means")
        new_df <- rbind(new_df, as.data.frame(mm))
      }
    }
    
  } # for loop
  return(new_df)
}

test_spar_pdfs <- function(comb_df) {
  spar_vals <- seq(.15,1.55, by=.05)
  pdf(file="/Users/ncolaian/Documents/phanstiel_lab/data/pdf_out/bin_spar_an50.pdf", w=11, h=8)
  smooth_min <- min(max(comb_df$means[comb_df$model == "Loop&FL"]),
                    max(comb_df$means[comb_df$model == "TADs"]),
                    max(comb_df$means[comb_df$model == "Background"]))
  
  for( i in spar_vals) {
    t_spline <- smooth.spline(comb_df$distance[comb_df$model == "TADs" & comb_df$distance <= smooth_min], comb_df$means[comb_df$model == "TADs" & comb_df$distance <= smooth_min], spar=i)
    lf_spline <- smooth.spline(comb_df$distance[comb_df$model == "Loop&FL" & comb_df$distance <= smooth_min], comb_df$means[comb_df$model == "Loop&FL" & comb_df$distance <= smooth_min], spar=i)
    b_spline <- smooth.spline(comb_df$distance[comb_df$model == "Background" & comb_df$distance <= smooth_min], comb_df$means[comb_df$model == "Background" & comb_df$distance <= smooth_min], spar=i)
    
    t_df <- data.frame(t_spline$x, t_spline$y)
    lf_df <- data.frame(lf_spline$x, lf_spline$y)
    b_df <- data.frame(b_spline$x, b_spline$y)
    colnames(t_df) <- c("distance", "means")
    colnames(lf_df) <- c("distance", "means")
    colnames(b_df) <- c("distance", "means")
    
    t_new_ys <- predict(t_spline, seq(1,200))
    t_df <- rbind(t_df, as.data.frame(t_new_ys,col.names = c("distance","means")))
    
    lf_new_ys <- predict(lf_spline, seq(1,200))
    lf_df <- rbind(lf_df, as.data.frame(lf_new_ys,col.names = c("distance","means")))
    
    b_new_ys <- predict(b_spline, seq(1,200))
    b_df <- rbind(b_df, as.data.frame(b_new_ys,col.names = c("distance","means")))
    
    comb2 <- rbind(cbind(t_df, model = "TADs"),cbind(lf_df, model = "Loop&FL"),cbind(b_df, model = "Background" ))
    
    par(mfrow=c(2,1))
    
    g1 <- ggplot( comb_df, aes( x = distance, y = means, col=model))+
      geom_line()+
      ggtitle("Distance Vs Mean - Real")+
      theme(plot.title = element_text(hjust = .5))
    
    g2 <- ggplot( comb2, aes( x = distance, y = means, col=model))+
      geom_line()+
      ggtitle(paste("Smooth Spline Spar =", i, "Distance Vs Mean"))+
      theme(plot.title = element_text(hjust = .5))
    
    print(g1)
    print(g2)
  }
  dev.off()
}

trial_pdfs <- function(comb_df) {
#want to find change of slope in background and apply that to the end of the tad and loops
  tmax <- max(comb_df$distance[ comb_df$model == "TADs" ])
  
}

#### MAIN ####
#TEST#
#This is to test the binned data anaysis on all the reads we could find associated with tads and loop/flares
#tad <- ddply(data4tad,"distance",summarize,total = sum(reads, na.rm = TRUE), means = mean(reads, na.rm = TRUE), sd = sd(reads, na.rm = TRUE), N = sum(!is.na(reads)), model = "TADs", tot_mean = total/N)
#fandl <- ddply(data4loop,"distance",summarize, total = sum(reads, na.rm = TRUE), means = mean(reads, na.rm = TRUE), sd = sd(reads, na.rm = TRUE), N = sum(!is.na(reads)), model = "Loop&FL", tot_mean = total/N)
#back <- ddply(data4back,"distance",summarize, total = sum(reads, na.rm = TRUE), means = mean(reads, na.rm = TRUE), sd = sd(reads, na.rm = TRUE), N = sum(!is.na(reads)), model = "Background", tot_mean = total/N)

#separate df
tad <- comb[comb$model=="TADs",]
fandl <- comb[comb$model=="Loop&FL",]
back <- comb[comb$model=="Background",]

#create binned df
tad <- combine_bins_if_necessary(450, tad )
fandl <- combine_bins_if_necessary(450, fandl)
back <- combine_bins_if_necessary(450, back)

#somehow got to matrix of lists? but bring it back to reg matrix
tad <- as.data.frame(sapply(tad, function(x){as.numeric(x)}))
fandl <- as.data.frame(sapply(fandl, function(x){as.numeric(x)}))
back <- as.data.frame(sapply(back, function(x){as.numeric(x)}))

#bring the dfs back together
tad$model <- "TADs"
fandl$model <- "Loop&FL"
back$model <- "Background"

new_comb <- rbind(tad,fandl,back)

test_spar_pdfs(new_comb)
trial_pdfs(new_comb)

