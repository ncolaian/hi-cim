#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(Sushi,lib.loc="/nas02/home/n/c/ncolaian/R/x86_64-pc-linux-gnu-library/3.2")
library(Sushi)# for testing
library(MASS)

#Handles input and puts the options passed into the program into a matrix
#Matrix makes it easier to look up values.
params = matrix(c(
  "loop_file", "l", 1, "character",
  "distr_data", "d", 1, "character",
  "chrom","c", 1, "integer",
  "start", "st", 1, "integer",
  "end", "e", 1, "integer",
  "bin_length", "b",1, "integer",
  "tad_distr", "t",1,"character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

#test
opt <- c("l", "a", 1, 1, 1, 1, "t","n", "norm_factor")
names(opt) <- c("loop_file","distr_data", "chrom", "start", "end", "bin_length", "tad_distr")
opt$loop_file <- "/Users/ncolaian/Documents/phanstiel_lab/data/CI_THP1_O_0.0.0.loops.10Kb.bedpe"
opt$distr_data <- "/Users/ncolaian/Documents/phanstiel_lab/data/distr_192.txt"
opt$chrom <- 20
opt$start <- 48500000
opt$end <- 50100000
opt$bin_length <- 10000
opt$tad_distr <- "/Users/ncolaian/Documents/phanstiel_lab/data/tad_distr_146.txt"
opt$tad_loop_dist <- "/Users/ncolaian/Documents/phanstiel_lab/data/tad_sig_loop_dist.txt"
opt$noise <- "/Users/ncolaian/Documents/phanstiel_lab/data/noise.txt"
opt$norm_factors <- "/Users/ncolaian/Documents/phanstiel_lab/data/CI_THP1_O_0.0.0.chr20.10Kb.MQ30.KR_norm.vals"



#### SUBROUTINES ####

#create a blank matrix that will hold the simulated data
create_sim_matrix <- function(starting, ending, bl) {
  options(scipen = 999)
  vec <- seq(starting, ending, by = bl)
  new_matrix <- matrix(0L, nrow = length(vec), ncol = length(vec))
  colnames(new_matrix) <- paste(vec,"lab", sep="_")
  rownames(new_matrix) <- paste(vec,"lab", sep="_")
  return(new_matrix)
}

get_loop_data <- function(starting, ending, loop, chr) {
  needed_vals <- loop[loop$bbin_start <= ending & loop$abin_start >= starting & as.character(loop$achr) == paste("chr",chr,sep=""),] #  
  return(needed_vals)
}

get_areas_in_tads_flares_loops <- function(need_loop, real, bl,flares=TRUE, tads=TRUE) {
  loop_vec <- c()
  flare_vec <- c()
  tad_vec <- c()
  tad_loop <- c()
  flare_loop <- c()
  
  for(i in seq(1:length(need_loop$abin_start))) {
    starting <- need_loop$abin_start[i]
    ending <- need_loop$bbin_start[i]
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

add_mean_value <- function(x,y, full_matrix, distr, tad_distr, tad_loop_dist, vec_list, bl, w_bg=TRUE) {
  options(scipen = 10, digits=10)
  num_loops <- sum((vec_list[[1]] == paste(x,":", y, sep = "")), na.rm = TRUE )
  if ( num_loops == 0 ){
    num_loops <- num_loops + sum((vec_list[[1]] == paste(y,":", x, sep = "")), na.rm = TRUE )
  }
  
  num_flares <- sum((vec_list[[2]] == paste(x,":", y, sep = "")), na.rm = TRUE)
  if ( num_flares == 0 ){
    num_flares <- num_flares + sum((vec_list[[2]] == paste(y,":", x, sep = "")), na.rm = TRUE)
  }
  
  num_tads <- sum((vec_list[[3]] == paste(x,":", y, sep = "")), na.rm = TRUE )
  if ( num_tads == 0 ){
    num_tads <- num_tads + sum((vec_list[[3]] == paste(y,":", x, sep = "")), na.rm = TRUE )
  }
  
  dist <- abs((y-x)/bl)
  if(dist>192){
    dist <- 192
  }
  
  value <- 0;
  if(num_loops > 0) {
    loop_val <- distr$mu[distr$distance == 0 & distr$model == "Loop&FL"]
    value <- value + loop_val
  }
  
  if(num_flares > 0) {
    bg_val <- distr$mu[distr$distance == dist & distr$model == "Background"]
    flare_val <- distr$mu[distr$distance == dist & distr$model == "Loop&FL"]
    factor_v <- (flare_val/bg_val)^num_flares
    
    value <- value + bg_val*factor_v
  }
  
  if(num_tads > 0) {
    #get the location in the tad_sig values are held
    tad_loops <- which(vec_list[[3]] %in% paste(x,":", y, sep = ""))
    tad_loops <- append(tad_loops, which(vec_list[[3]] %in% paste(y,":", x, sep = "")))
    
    #get the min loop pixel is associated with
    loop <- min(vec_list[[5]][tad_loops])
    if (dist > max(tad_distr$distance) ) {
      dist <- max(tad_distr$distance)
    }
    t_val <- 0
    if(num_tads == 1) {
      sig_fact <- tad_signal_vs_loop(loop)
      t_val <- t_val + sig_fact * distr$mu[distr$distance == dist & distr$model == "TADs"]
    }
    else if ( num_tads == 2) {
      #Take into account of amount of loops
      t_fact <- (tad_distr$mu[tad_distr$distance == dist & tad_distr$model == "TAD2"]/
                   tad_distr$mu[tad_distr$distance == dist & tad_distr$model == "TADs"])
      sig_fact <- tad_signal_vs_loop(loop)
      #take into account loop dist associated with pixel
      t_val <- t_val + t_fact*sig_fact*distr$mu[distr$distance == dist & distr$model == "TADs"]
    }
    else {
      #Take into account of amount of loops
      t_fact <- (tad_distr$mu[tad_distr$distance == dist & tad_distr$model == "TAD_more"]/
                   tad_distr$mu[tad_distr$distance == dist & tad_distr$model == "TADs"])
      #take into account loop dist associated with pixel
      if(t_fact < 1) {
        t_fact <- 1/t_fact
      }
      sig_fact <- tad_signal_vs_loop(loop)
      t_val <- t_val + t_fact*sig_fact*distr$mu[distr$distance == dist & distr$model == "TADs"]
    }
    
    value <- value + t_val
    dist <- abs((y-x)/bl)
  }
  
  #this means it is a background pixel
  if(value == 0 && w_bg == TRUE) {
    val <- distr$mu[distr$distance == dist & distr$model == "Background"]
    value <- value + val
  }
  #val_vec <- rgamma(1, as.numeric(distr$scale[distr$distance == dist & distr$model == "Background"]), rate = as.numeric(distr$rate[distr$distance == dist & distr$model == "Background"]))
  #value <- value + sum(as.integer(val_vec))
  
  #test for the drawing of only one feature
  #value <- value + max(as.integer(val_vec))
  
  if(is.na(value)){
    value <- 0
  }
  
  return(value)
}

tad_signal_vs_loop <- function(x) {
  val <- .258379 + .176503*log(x) + rnorm(1, mean = 0, sd = .1169^2)
  return(val)
}

create_hic_heatmap <- function(matrix1, matrix2, chr, starting, ending, bl, zr=c(0,400),
                               legend=T, lab="Simulated Data" ) {
  chromo <- paste("chr", chr, sep="")
  vec <- seq(starting, ending, by = bl)
  colnames(matrix1) <- vec
  rownames(matrix1) <- vec
  colnames(matrix2) <- vec
  rownames(matrix2) <- vec
  
  hcplot <- plotHic2(matrix1, chromo, starting, ending, zrange=zr, palette = colorRampPalette(c("white","dodgerblue2")), labeltype = "bp", resolution = 10000, plottype = "square", half = "top", format = "full")
  hcplot <- plotHic2(matrix2, chromo, starting, ending, zrange=zr, palette = colorRampPalette(c("white","firebrick2")), labeltype = "bp", resolution = 10000, plottype = "square", half = "bottom", format = "full", add=T)
  labelgenome(chrom=chromo,starting, ending, scale="Mb")
  labelgenome(chrom=chromo,starting, ending, scale="Mb", side = 2, las = 2)
  title(main="Simulated Hi-C Data")
  
  return(hcplot)
}

compare_real_hic_heatmap <- function(matrix1, real_mat, chr, starting, ending, bl, zr=c(0,300),
                                     legend=T, lab="Simulated Data", sim_v_real=TRUE ) {
  chromo <- paste("chr", chr, sep="")
  
  if ( sim_v_real == FALSE ) {
    hcplot <- plotHic2(matrix1, chromo, starting, ending, zrange=zr, palette = colorRampPalette(c("white","firebrick2")), labeltype = "bp", resolution = 10000, plottype = "square", half = "top", format = "sparse")
  }
  else {
    vec <- seq(starting, ending, by = bl)
    colnames(matrix1) <- vec
    rownames(matrix1) <- vec
    hcplot <- plotHic2(matrix1, chromo, starting, ending, zrange=zr, palette = colorRampPalette(c("white","firebrick2")), labeltype = "bp", resolution = 10000, plottype = "square", half = "top", format = "full")
  }
  hcplot <- plotHic2(real_mat, chromo, starting, ending, zrange=zr, palette = colorRampPalette(c("white","firebrick2")), labeltype = "bp", resolution = 10000, plottype = "square", half = "bottom", format = "sparse", add=T)
  labelgenome(chrom=chromo,starting, ending, scale="Mb")
  labelgenome(chrom=chromo,starting, ending, scale="Mb", side = 2, las = 2)
  title(main="Real Hi-C Data")
  
  return(hcplot)
}

denormalize_sample_and_save_new_val <- function(x,y,value,norm_factor,noise_fxn,bl) {
  #Need to add the +1 since the bins start from 0 but indexing starts at 1
  if(!(x == 0)) {
    x <- x/bl + 1
  }
  else {
    x <- x+1
  }
  y <- y/bl + 1
  
  #Denormalize value  
  value <- denormalize(value, x, y, norm_factor)
  
  #Get value considering signal noise
  value <- rnorm(1, mean = value, sd = sqrt(noise_fxn(value)) )
  
  #renormalize
  value <- normalize(value, x,y, norm_factor)
  
  return(value)
}

#norm factors is a vector of factors assocated with each bin
normalize <- function(prenorm, p1, p2, normFactors){
  
  # normalize based on normalization factors and distances
  prenorm.norm = prenorm/(normFactors[p1,]*normFactors[p2,])
  
  return(prenorm.norm)
}

denormalize <- function(norm, p1, p2, normFactors){
  
  # normalize based on normalization factors and distances
  norm.denorm = normFactors[p1,]*normFactors[p2,]*norm
  
  return(norm.denorm)
}


#### MAIN ####
#loop data
loop <- read.delim(opt$loop_file, header = FALSE)
loop <- as.data.frame(loop)
colnames(loop) <- c("achr", "abin_start", "abin_end", "bchr", "bbin_start", "bbin_end",
                    "bin_length", "loop_pval")
#assign new bin length values to make it a number and not a character
loop$distance <- abs(loop$abin_start - loop$bbin_start)

#distribution data
distr_data <- read.delim(opt$distr_data)
distr_data <- as.data.frame(distr_data)
distr_data$model <- as.character(distr_data$model)

#tad_dist
tad_distr_data <- read.delim(opt$tad_distr)
tad_distr_data$model <- as.character(tad_distr_data$model)

#tad_loop_dist
tad_loop_dist <- read.delim(opt$tad_loop_dist)
tad_loop_dist$model <- as.integer(tad_loop_dist$model)

#noise_file
noise <- read.delim(opt$noise)
noise_fxn <- splinefun(noise)

#norm_factors
nfact <- read.delim(opt$norm_factors)

#real_data
real_data <- read.delim("/Users/ncolaian/Documents/phanstiel_lab/data/count_matrix_chr20_norm.txt", header = TRUE)


#bin_length fix
#opt$bin_length <- opt$bin_length * 1000

#sim_mat <- create_sim_matrix(opt$start, opt$end, opt$bin_length)

loops_to_sim <- get_loop_data(opt$start, opt$end, loop, opt$chrom)

non_bg <- get_areas_in_tads_flares_loops(loops_to_sim, real_data, opt$bin_length, flares = FALSE, tads = TRUE)

#go through and add values to the sim matrix
#for(i in rownames(sim_mat)) {
#  sim_mat[i,] <- sapply(colnames(sim_mat), function(x) { add_value(as.integer(gsub( "_lab",replacement = "", x)), as.integer(gsub("_lab",replacement = "",i)), sim_mat, distr_data, non_bg, opt$bin_length)} )
#}

# mean values used
mean_mat <- create_sim_matrix(opt$start, opt$end, opt$bin_length)
for(i in rownames(mean_mat)) {
  mean_mat[i,] <- sapply(colnames(mean_mat), function(x) { add_mean_value(as.integer(gsub( "_lab",replacement = "", x)), as.integer(gsub("_lab",replacement = "",i)), mean_mat, distr_data, tad_distr_data,tad_loop_dist, non_bg, opt$bin_length)} )
}

#Normalize Data
for(i in rownames(mean_mat)) {
  mean_mat[i,] <- sapply(colnames(mean_mat), function(x) { denormalize_sample_and_save_new_val(as.integer(gsub( "_lab",replacement = "", i)), as.integer(gsub("_lab",replacement = "",x)), mean_mat[i,x],nfact, noise_fxn, opt$bin_length)} )
}

#Get the actual data
#opt$start <- 48700000
#opt$end <- 49200000
figure_mat <- real_data[real_data$Bin1 >= opt$start & real_data$Bin2 <= opt$end, c("Bin1", "Bin2", "Signal") ]

#hicp <- create_hic_heatmap(mean_mat,sim_mat, opt$chrom, opt$start, opt$end, opt$bin_length)
compare_hicp <- compare_real_hic_heatmap(mean_mat,figure_mat, opt$chrom, opt$start, opt$end, opt$bin_length)
compare_hicp <- compare_real_hic_heatmap(figure_mat,figure_mat, opt$chrom, opt$start, opt$end, opt$bin_length, sim_v_real = FALSE)

png("sim_hic.png")
print(hicp)
dev.off()


