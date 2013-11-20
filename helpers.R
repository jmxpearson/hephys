# helper functions for R analysis
source('setup_env.R')

# function to calculate roc curve
roc <- function(fitobj, y, X, dth=0.01) {
  # make a dataframe of true positive, false positive measures for a
  # variety of thresholds
  thvec <- seq(0, 1, dth)
  N <- length(thvec)
  TPR <- rep(NA, N)
  FPR <- rep(NA, N)
  for (ind in 1:N) {
    pp <- predict(fitobj, X, type='response') > thvec[ind]
    TP <- sum(pp == TRUE & y == TRUE)
    FP <- sum(pp == TRUE & y == FALSE)
    TN <- sum(pp == FALSE & y == FALSE)
    FN <- sum(pp == FALSE & y == TRUE)
    TPR[ind] <- TP/(TP + FN)
    FPR[ind] <- FP/(FP + TN)
  }
  roc <- data.frame(thresh=thvec, TPR=TPR, FPR=FPR)
  return (roc)
}

# function to return dataframe suitable for plotting roc curve
get_performance <- function(ind) {
    glmobj <- fitobjs[[ind]]$glmobj
    fname <- paste(paste(chanlist[ind,], collapse='.'), 'lfpglmdata.csv', sep='.')
    dfile <- paste(ddir, fname, sep='/')
    dat <- read.csv(dfile, header=TRUE, row.names=1)

    # get performance
    min.ind <- which(glmobj$lambda == glmobj$lambda.min)
    auc <- glmobj$cvm[min.ind]

    # make predictions, calculate ROC
    y <- dat[,1]
    X <- as.matrix(dat[,-1])
    perf <- roc(glmobj, y, X, dth=0.001)
    return (list(perf=perf, auc=auc))
} 

# plot an roc curve
plotroc <- function(perf) {
    titlestr = paste('Area Under the Curve =', sprintf('%.2f', perf$auc))
    
    plt <- ggplot(perf$perf, aes(FPR, TPR))
    plt + geom_line(size=1.5) + geom_abline(slope=1, intercept=0) +
      labs(title = titlestr, x='\nFalse Positive Rate', y='True Positive Rate\n') +
      scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
      scale_y_continuous(limits=c(0,1), expand=c(0,0)) + theme_bw() +
      theme(plot.title = element_text(lineheight=1, size=36), 
            axis.text=element_text(color='black', size=12),
            axis.title.x=element_text(size=20), axis.title.y=element_text(size=20)
            )
}