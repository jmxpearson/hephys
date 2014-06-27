# set up file structure and common variables used by R code

adir <- '~/code/hephys/bartc'  # analysis directory
ddir <- '~/data/bartc'  # data directory

# lfp data
cfile <- 'data/lfp_channel_file.csv'
chanlist <- read.csv(paste(adir, cfile, sep='/'), header=FALSE)
chanlist <- chanlist[, 1:2]
chanlist <- chanlist[!duplicated(chanlist),]
numlfp <- dim(chanlist)[1]

# spike data
cfile <- 'data/valid_units.csv'
unitlist <- read.csv(paste(adir, cfile, sep='/'), header=FALSE)
numunits <- dim(unitlist)[1]

library(glmnet)
library(ggplot2)
library(plyr)
library(reshape)