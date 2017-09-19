#! /usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

#the combined data points
data4distr <- read.delim(args[1])

