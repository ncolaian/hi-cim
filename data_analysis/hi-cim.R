#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(Sushi)

#Handles input and puts the options passed into the program into a matrix
#Matrix makes it easier to look up values.
params = matrix(c(
  "loop_file", "l", 1, "character",
  "distr_data", "d", 1, "character",
  "chrom","c", 1, "integer",
  "start", "s", 1, "integer",
  "end", "e", 1, "integer",
  "bin_length", "b", "integer"
), byrow = TRUE, ncol = 4)
opt = getopt(params)


#### SUBROUTINES ####

#create a blank matrix that will hold the simulated data
create_sim_matrix <- function(starting, ending, bl) {
  vec <- seq(starting, ending, by = bl)
  new_matrix <- matrix(0L, nrow = length(vec), ncol = length(vec))
  colnames(new_matrix) <- vec
  rownames(new_matrix) <- vec
  return(new_matrix)
}

get_loop_data <- function(starting, ending, loop, chr) {
  needed_vals <- loop[loop$abin_start >= starting & loop$bbin_start <= ending & loop$achr == chr,]
  return(needed_vals)
}

get_areas_in_tads_flares_loops <- function(need_loop, bl) {
  loop_vec <- c()
  flare_vec <- c()
  tad_vec <- c()

  for(i in length(need_loop)) {
    starting <- need_loop$abin_start[i]
    ending <- need_loop$bbin_start[i]
    distance <- ending - starting
    
    loop_vec <- c(loop_vec, paste(starting, ":", ending, sep = ""))
    flare_vec <- c(flare_vec, get_flares(starting, ending, bl))
    tad_vec <- c(tad_vec, get_tads(starting,ending, bl))
  }
  
  return(list(loop_vec, flare_vec, tad_vec))
}

get_flares <- function(starting, ending, bl) {
  vec<- c()
  adjust <- bl*2
  flare_end <- endind+adjust
  bot <- seq(flare_start, starting+adjust, by = bl)
  top <- seq(flare_end, ending-adjust, by = bl)
  
  
  #get left flare and top right combined
  for(spot in bot) {
    new <- vapply(seq(spot, flare_end, by = bl), function(x) {paste(spot,":",x, sep="")})
    vec <- c(vec, new)
  }
  
  #get top flare
  for(spot in top) {
    new <- vapply(seq(spot, (starting-adjust), by = bl), function(x) {paste(spot,":",x, sep = "")})
    vec <- c(vec, new)
  }
  
  return(vec)
}

get_tads <- function(starting, ending, bl) {
  vec <- c()
  adjust <- bl*2
  new_s <- starting+adjust
  new_e <- ending-adjust
  
  for(spot in seq(new_s, new_e, by=bl)) {
    new <- vapply(seq(spot, new_e, by = bl), function(x) {paste(spot, ":", x, sep = "")})
    vec <- c(vec, new)
  }
  
  return(vec)
}

add_value <- function(x,y, full_matrix, distr, vec_list, bl) {
  num_loops <- length(vec_list[[1]] == paste(x,":", y, sep = ""))
  num_flares <- length(vec_list[[2]] == paste(x,":", y, sep = ""))
  num_tads <- length(vec_list[[3]] == paste(x,":", y, sep = ""))
  
  dist <- (y-x)/bl
  value <- 0;
  if(num_loops != 0) {
    val_vec <- rgamma(num_loops, distr$shape[distr$distance == 1, distr$model == "Loop&Fl"], rate = distr$rate[distr$distance == 1, distr$model == "Loop&Fl"])
    value <- value + sum(val_vec)
  }
  
  if(num_flares != 0) {
    val_vec <- rgamma(num_flares, distr$shape[distr$distance == dist, distr$model == "Loop&Fl"], rate = distr$rate[distr$distance == dist, distr$model == "Loop&Fl"])
    value <- value + sum(val_vec)
  }
  
  if(num_tads != 0) {
    val_vec <- rgamma(num_tads, distr$shape[distr$distance == dist, distr$model == "TADs"], rate = distr$rate[distr$distance == dist, distr$model == "TADs"])
    value <- value + sum(val_vec)
  }
  
  #this means it is a background pixel
  if(value == 0) {
    val_vec <- rgamma(1, distr$shape[distr$distance == dist, distr$model == "Background"], rate = distr$rate[distr$distance == dist, distr$model == "Backgroud"])
    value <- sum(val_vec)
  }
  
  full_matrix[as.character(x), as.character(y)] <- value
}

create_hic_heatmap <- function(matrix, chr, starting, ending, zr=c(0,200), pal=pal_O,
                               legend=T, lab="Simulated Data", res=hicres) {
  chromo <- paste("chr", chr, sep="")
  
  plotHic(matrix, chromo, starting, ending, zrange=zr, palette = pal, resolution=res)
  labelgenome(chrom=chromo,starting, ending, scale="Mb")
  title(main="Simulated Hi-C Data")
}


#### MAIN ####
#loop data
loop <- read.delim(opt$loop_file, header = FALSE)
loop_df <- as.data.frame(loop)
colnames(loop_df) <- c("achr", "abin_start", "abin_end", "bchr", "bbin_start", "bbin_end",
                       "bin_length", "loop_pval")
#assign new bin length values to make it a number and not a character
loop_df$distance <- abs(loop_df$abin_start - loop_df$bbin_start)

#norm count data
read_counts <- read.delim(opt$count_file, header = FALSE)
read_counts <- as.data.frame(read_counts)
colnames(read_counts) <- c("start", "end", "reads")

#distribution data
distr_data <- read.delim(opt$distr_data)
distr_data <- as.data.frame(distr_data)

#bin_length fix
opts$bin_length <- opts$bin_length * 10000




