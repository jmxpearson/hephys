fit_all_and_save <- function(filext, outname, family, datalist, measure, lambdatype) {
  fitobjs <- list()
  for (ind in 1:dim(datalist)[1]) {
    fname <- paste(paste(datalist[ind,], collapse='.'), filext, sep='.')
    dfile <- paste(ddir, fname, sep='/')
    print(dfile)
    dat <- read.table(dfile, sep=',', header=TRUE, row.names=1, colClasses=c('numeric'))
    thisfit <- run_glm(dat, family, measure, lambdatype)
    fitobjs[[ind]] <- thisfit
  }

  save(fitobjs, file=paste(ddir, outname, sep='/'))
}

run_glm <- function(dframe, type='binomial', measure="deviance", lambdatype='1se'){
  # given an input data frame, perform an elastic net regression
  # on the data contained therein

  # make model matrix
  y <- dframe[, 1]
  X <- dframe[, -1]

  allobjs = list()
  alphalist = seq(0, 1, 0.1)

  folds <- partition_data(y, nfolds=10)
  
  # loop over alpha values
  for (alpha in alphalist) {
    # run sparse regression
    # don't include constant in X; glmnet will fit an intercept, 
    # but if we supply
    # one manually, it will be penalized, which we don't want
    glmobj <- cv.glmnet(as.matrix(X), y, alpha = alpha, family = type, 
                        foldid = folds, intercept = TRUE, type.measure = measure)

    allobjs[[length(allobjs) + 1]] <- get_best_beta(glmobj, lambdatype)
  }

  bestobj <- get_best_alpha(allobjs, alphalist)

  plot(bestobj$glmobj)
  fit <- bestobj$glmobj$glmnet.fit
  plot(fit, xvar='lambda', label=TRUE)

  return(bestobj)
}

partition_data <- function(outcome, nfolds) {
  suppressMessages(require('caret'))

  folds <- createFolds(outcome, nfolds)

  foldid <- rep(NA, length(outcome))

  for (ind in 1:nfolds) {
    foldid[folds[[ind]]] <- ind
  }

  return(foldid)
}

get_best_beta <- function(glmobj, lambdatype) {
  # given a glmobj from elastic net, extract the cross-validation score and
  # coefficients of the best fit

  if (lambdatype == '1se') { minlambda <- glmobj$lambda.1se }
  else if (lambdatype == 'min') { minlambda <- glmobj$lambda.min }

  minlambda.ind <- which(glmobj$lambda == minlambda)
  fit <- glmobj$glmnet.fit
  beta <- fit$beta[, minlambda.ind]
  intercept <- fit$a0[minlambda.ind]
  score <- glmobj$cvm[minlambda.ind]

  return(list(beta = beta, intercept = intercept, score = score, glmobj = glmobj))
}

get_best_alpha <- function(objlist, alphalist) {
  # given a list of objects returned by get_best_beta and a list of alpha 
  # values, return the object corresponding to the highest score

  measure_type <- tolower(objlist[[1]]$glmobj$name)
  if (grepl('auc', measure_type)) {
    extract_best <- which.max
  } else {
    extract_best <- which.min
  }

  scorelist <- sapply(objlist, function(x) {x$score})
  best_ind <- extract_best(scorelist)  # if using auc, want max score

  bestobj <- c(objlist[[best_ind]], alpha=alphalist[best_ind])

  return(bestobj)
}