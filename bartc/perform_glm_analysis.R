suppressMessages(library(glmnet))
source('glm_helpers.R')
source('setup_env.R')

args <- commandArgs(TRUE)

if (args[1] == 'lfp') {
    print("Performing LFP analysis...")
    filext <- 'lfpglmdata.csv'
    outname <- 'lfpfitdata' 
    family <- 'binomial'
    datalist <- chanlist
    measure <- 'auc'
} else if (args[1] == 'spikes') {
    print("Performing spike analysis...")
    filext <- 'spkglmdata.csv'
    outname <- 'spkfitdata'
    family <- 'poisson'
    datalist <- unitlist
    measure <- 'deviance'
}

fit_all_and_save(filext, outname, family, datalist, measure)