library(glmnet)
source('dolfpglm.R')

adir <- '~/code/hephys'  # analysis directory
ddir <- '~/data/bartc'  # data directory
cfile <- 'lfp_channel_file.csv'
chanlist <- read.csv(paste(adir, cfile, sep='/'), header=FALSE)
chanlist <- chanlist[, 1:2]
chanlist <- chanlist[!duplicated(chanlist),]
numunits <- dim(chanlist)[1]

fitobjs = list()
for (ind in 1:numunits) {
  fname <- paste(paste(chanlist[ind,], collapse='.'), 'lfpglmdata.csv', sep='.')
  dfile <- paste(ddir, fname, sep='/')
  print(dfile)
  thisfit <- dolfpglm(dfile)
  fitobjs[[ind]] <- thisfit
}

save(fitobjs, file=paste(ddir, 'lfpfitdata', sep='/'))