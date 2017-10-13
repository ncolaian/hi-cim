#This script will attempt to classify tads as strong vs weak tads
#Will also save the loop signal

library(MASS)
install.packages("/Users/ncolaian/Library/R/3.3/libraryphanstielR")
library(phanstielR)
library(ggplot2)

#### SUBROUTINES ####
get_indiv_tad_matrix <- function(pix_row, full_pix, tad_dist) {
  full_pix <- full_pix[full_pix$Bin1 >= pix_row$Bin1[1] & full_pix$Bin2 <= pix_row$Bin2[1] &
                         !(full_pix$Bin1 == pix_row$Bin1[1] && full_pix$Bin2 == pix_row$Bin2[1]),]
  
  full_pix$Dist <- abs(full_pix$Bin2 - full_pix$Bin1)/10000
  #need to normalize the signal
  for(i in 1:nrow(full_pix)) {
    full_pix$Signal[i] <- full_pix$Signal[i]/tad_dist$mu[tad_dist$distance == full_pix$Dist[i]]
  }
  full_pix <- full_pix[,c("Dist", "Signal")]
  return(full_pix)
}

get_distributed_tads <- function(tad_vals, tad_dist) {
  dist_vs_counts_tads <- matrix(ncol = 2)
  dist_vs_counts_tads <-as.data.frame(dist_vs_counts_tads)
  colnames(dist_vs_counts_tads) <- c("Dist", "Signal")
  
  for(i in seq(0,max(tad_vals$Dist))) {
    if(as.integer(nrow(tad_vals[tad_vals$Dist == i,]) <= 1)) next
    val_vec <- rgamma(nrow(tad_vals[tad_vals$Dist == i,]),tad_dist$scale[tad_dist$distance == i], rate = tad_dist$rate[tad_dist$distance == i])
    val_vec <- val_vec/tad_dist$mu[tad_dist$distance == i]
    mat <- data.frame(rep(i, length(val_vec)), val_vec)
    dist_vs_counts_tads <- rbind.phan(list(dist_vs_counts_tads,mat))
  }
  
  return(dist_vs_counts_tads)
}

#### MAIN ####
distr <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/distr_192.txt")
pixel_vals <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/count_matrix_chr20.txt")

#Get the loop coordinates and create a list 
loop_pix <- pixel_vals[pixel_vals$Loops > 0,]

#get tad distr
t_distr <- distr[distr$model == "TADs",]

#declare a list that will contain the pixel values for each list
tad_vals <- vector("list", length = nrow(loop_pix))

#declare a list for the sampled vals
samp_vals <- vector("list", length = nrow(loop_pix))

#declare_t_test_list
t_test_l <- vector("list", length = nrow(loop_pix))

#go through and get tad matrices
for(i in seq(1,(nrow(loop_pix)))) {
  print(i)
  tad_vals[[i]] <- get_indiv_tad_matrix(loop_pix[i,], na.omit(pixel_vals), t_distr)
  samp_vals[[i]] <- get_distributed_tads(tad_vals[[i]], t_distr)
  t_test_l[[i]] <- t.test(tad_vals[[i]], samp_vals[[i]])
}
vec <- c()
for(i in seq(1,(nrow(loop_pix)))) {
  if(t_test_l[[i]]$p.value < .5) {
    print(i)
    vec <- append(vec, i)
    print(t_test_l[[i]]$p.value)
  }
}
str(samp_vals[[1]])
signal_vec_low <- c()
signal_vec_high <- c()
for( i in vec ) {
  if ( median(tad_vals[[i]]$Signal) < median(samp_vals[[i]]$Signal, na.rm = TRUE) ) {
    signal_vec_low <- append(signal_vec_low, loop_pix$Signal[i])
  }
  else {
    signal_vec_high <- append(signal_vec_high, loop_pix$Signal[i])
  }
}
norm_sig <- c()
for ( i in seq(1:length(loop_pix$Signal)) ) {
  if ( !(i %in% vec) ) {
    norm_sig <- append(norm_sig, loop_pix$Signal[i])
  }
}
mm_high <- data.frame(signal_vec_high)
mm_norm <- data.frame(norm_sig)
mm_high$model <- "high"
mm_norm$model <- "norm"

diff_mat <- rbind.phan(list(mm_high, mm_norm))
colnames(diff_mat) <- c("Signal", "Model")


ggplot(diff_mat, aes(Model,Signal))+
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter()+
  ggtitle(label = "Loop Signal and Tad Strength")+
  theme( plot.title = element_text(hjust = .5))

t.test(diff_mat$Signal[diff_mat$Model == "norm"], diff_mat$Signal[diff_mat$Model == "high"])

#Plot mean tad strength vs loop strength
tad_sig <- c()
for ( i in seq(1:length(loop_pix$Signal)) ) {
  tad_sig <- append(tad_sig, median(tad_vals[[i]]$Signal))
}
tad_v_loop <- data.frame(tad_sig, loop_pix$Signal, ((loop_pix$Bin2-loop_pix$Bin1)/10000))
colnames(tad_v_loop) <- c("Tad_Signal", "Loop_Signal", "Distance")

ggplot(tad_v_loop, aes(Tad_Signal, log(Loop_Signal/Distance)))+
  geom_line()

cor(tad_v_loop$Tad_Signal,y = tad_v_loop$Loop_Signal/tad_v_loop$Distance, method = "pearson")
tvl <- lm(Tad_Signal ~ log(Loop_Signal/Distance), tad_v_loop)
summary(tvl)
par(mfrow=c(2,2))
plot(tvl) 
plot(tvl$fitted.values, tad_v_loop$Loop_Signal)
hist(tvl$residuals, breaks = 50)
dev.off()

#look at pixel value to distance
dist_v_pix <- data.frame(loop_pix$Signal,((loop_pix$Bin2-loop_pix$Bin1)/10000))
colnames(dist_v_pix) <- c("Loop_Signal", "Distance")

ggplot(dist_v_pix, aes(Distance, Loop_Signal))+
  geom_line()
cor(dist_v_pix, method = "pearson")
dvl <- lm(Loop_Signal ~ Distance, dist_v_pix)
plot(dvl)
