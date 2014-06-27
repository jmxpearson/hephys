library(glmnet)
source('glm_helpers.R')
source('setup_env.R')


fitobjs = list()
for (ind in 1:numspk) {
  fname <- paste(paste(chanlist[ind,], collapse='.'), 'spkglmdata.csv', sep='.')
  dfile <- paste(ddir, fname, sep='/')
  print(dfile)
  dat <- read.table(dfile, sep=',', header=TRUE, row.names=1, colClasses=c('numeric'))
  thisfit <- run_glm(dat, 'poisson')
  fitobjs[[ind]] <- thisfit
}

save(fitobjs, file=paste(ddir, 'spkfitdata', sep='/'))