#!/usr/bin/env Rscript

library(getopt)


#Handles input and puts the options passed into the program into a matrix
#Matrix makes it easier to look up values.
params = matrix(c(
  "loop_file", "f", 1, "character",
  "alt_loops", "f", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

