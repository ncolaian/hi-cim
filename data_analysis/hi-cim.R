#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(Sushi,lib.loc="/nas02/home/n/c/ncolaian/R/x86_64-pc-linux-gnu-library/3.2")
library(Sushi)# for testing

#Handles input and puts the options passed into the program into a matrix
#Matrix makes it easier to look up values.
params = matrix(c(
  "loop_file", "l", 1, "character",
  "distr_data", "d", 1, "character",
  "chrom","c", 1, "integer",
  "start", "st", 1, "integer",
  "end", "e", 1, "integer",
  "bin_length", "b",1, "integer"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

#test
opt <- c("l", "a", 1, 1, 1, 1)
names(opt) <- c("loop_file","distr_data", "chrom", "start", "end", "bin_length")
opt$loop_file <- "/Users/ncolaian/Documents/phanstiel_lab/data/CI_THP1_O_0.0.0.loops.10Kb.bedpe"
opt$distr_data <- "/Users/ncolaian/Documents/phanstiel_lab/data/distr_50.txt"
opt$chrom <- 20
opt$start <- 1600000
opt$end <- 4000000
opt$bin_length <- 10



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

get_areas_in_tads_flares_loops <- function(need_loop, bl) {
  loop_vec <- c()
  flare_vec <- c()
  tad_vec <- c()

  for(i in seq(1:length(need_loop$abin_start))) {
    starting <- need_loop$abin_start[i]
    ending <- need_loop$bbin_start[i]
    
    loop_vec <- append(loop_vec, paste(starting, ":", ending, sep = ""))
    flare_vec <- append(flare_vec, get_flares(starting, ending, bl))
    tad_vec <- append(tad_vec, get_tads(starting,ending, bl))
  }
  
  return(list(loop_vec, flare_vec, tad_vec))
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

get_tads <- function(starting, ending, bl) {
  vec <- c()
  adjust <- bl*2
  new_s <- starting+adjust
  new_e <- ending-adjust

  for(spot in seq(from = new_s, to = new_e, by=bl)) {
    new <- sapply(seq(spot, new_e, by = bl), function(x) {x <- paste(spot,":",x, sep=""); return(x)})
    vec <- append(vec, new)
  }
  
  return(vec)
}

add_value <- function(x,y, full_matrix, distr, vec_list, bl) {
  options(scipen = 10, digits=10)
  num_loops <- sum((vec_list[[1]] == paste(x,":", y, sep = "")), na.rm = TRUE )
  num_flares <- sum((vec_list[[2]] == paste(x,":", y, sep = "")), na.rm = TRUE)
  num_tads <- sum((vec_list[[3]] == paste(x,":", y, sep = "")), na.rm = TRUE )

  dist <- abs((y-x)/bl)
  if(dist>50){
    dist <- 50
  }
    
  value <- 0;
  if(num_loops > 0) {
    val_vec <- rgamma(num_loops, as.numeric(distr$scale[distr$distance == 0 & distr$model == "Loop&FL"]), rate = as.numeric(distr$rate[distr$distance == 0 & distr$model == "Loop&FL"]))
    value <- value + sum(round(val_vec))*2
  }
  
  if(num_flares > 0) {
    val_vec <- rgamma(num_flares, as.numeric(distr$scale[distr$distance == dist & distr$model == "Loop&FL"]), rate = as.numeric(distr$rate[distr$distance == dist & distr$model == "Loop&FL"]))
    value <- value + sum(round(val_vec))
  }
  
  if(num_tads > 0) {
    val_vec <- rgamma(num_tads, as.numeric(distr$scale[distr$distance == dist & distr$model == "TADs"]), rate = as.numeric(distr$rate[distr$distance == dist & distr$model == "TADs"]))
    value <- value + sum(round(val_vec))
  }
  
  #this means it is a background pixel
  # if(value == 0) {
  #   val_vec <- rgamma(1, as.numeric(distr$scale[distr$distance == dist & distr$model == "Background"]), rate = as.numeric(distr$rate[distr$distance == dist & distr$model == "Background"]))
  #   value <- sum(as.integer(val_vec))
  # }
  val_vec <- rgamma(1, as.numeric(distr$scale[distr$distance == dist & distr$model == "Background"]), rate = as.numeric(distr$rate[distr$distance == dist & distr$model == "Background"]))
  value <- value + sum(as.integer(val_vec))
  if(is.na(value)){
    value <- 0
  }
  
  return(value)
}

create_hic_heatmap <- function(matrix, chr, starting, ending, bl, zr=c(0,600),
                               legend=T, lab="Simulated Data" ) {
  chromo <- paste("chr", chr, sep="")
  vec <- seq(starting, ending, by = bl)
  colnames(matrix) <- vec
  rownames(matrix) <- vec
  
  hcplot <- plotHic2(matrix, chromo, starting, ending, zrange=zr, palette = colorRampPalette(c("white","dodgerblue2")), labeltype = "bp", resolution = 10000, plottype = "square", half = "both", format = "full")
  labelgenome(chrom=chromo,starting, ending, scale="Mb")
  labelgenome(chrom=chromo,starting, ending, scale="Mb", side = 2, las = 2)
  title(main="Simulated Hi-C Data")
  
  return(hcplot)
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

#bin_length fix
opt$bin_length <- opt$bin_length * 1000

sim_mat <- create_sim_matrix(opt$start, opt$end, opt$bin_length)

loops_to_sim <- get_loop_data(opt$start, opt$end, loop, opt$chrom)

non_bg <- get_areas_in_tads_flares_loops(loops_to_sim, opt$bin_length)

#go through and add values to the sim matrix
for(i in rownames(sim_mat)) {
  sim_mat[i,] <- sapply(colnames(sim_mat), function(x) { add_value(as.integer(gsub( "_lab",replacement = "", x)), as.integer(gsub("_lab",replacement = "",i)), sim_mat, distr_data, non_bg, opt$bin_length)} )
}

#new_sim_mat <- sapply( rownames(sim_mat), function(x) sapply(colnames(sim_mat), function(y) { 
    #add_value(as.integer(gsub("_lab",replacement = "",x)),as.integer(gsub("_lab",replacement = "" ,y)), sim_mat, distr_data, non_bg, opt$bin_length) }))

hicp <- create_hic_heatmap(sim_mat, opt$chrom, opt$start, opt$end, opt$bin_length)

png("sim_hic.png")
print(hicp)
dev.off()


