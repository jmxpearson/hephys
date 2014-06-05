library(glmnet)
source('dolfpglm.R')
source('setup_env.R')


fitobjs = list()
for (ind in 1:numunits) {
  fname <- paste(paste(chanlist[ind,], collapse='.'), 'lfpglmdata.csv', sep='.')
  dfile <- paste(ddir, fname, sep='/')
  print(dfile)
  thisfit <- dolfpglm(dfile)
  fitobjs[[ind]] <- thisfit
}

save(fitobjs, file=paste(ddir, 'lfpfitdata', sep='/'))