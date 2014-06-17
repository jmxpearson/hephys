
run_lfp_glm <- function(dfile){
  # given an input data file, perform an elastic net regression
  # on the data contained therein

  dat <- read.table(dfile, sep=',', header=TRUE, row.names=1, colClasses=c('numeric'))
  
  # make model matrix
  y <- dat[, 1]
  X <- dat[, -1]

  allobjs = list()
  alphalist = seq(0, 1, 0.1)

  folds <- partition_data(y, nfolds=10)
  
  # loop over alpha values
  for (alpha in alphalist) {
    # run sparse regression
    # don't include constant in X; glmnet will fit an intercept, 
    # but if we supply
    # one manually, it will be penalized, which we don't want
    glmobj <- cv.glmnet(as.matrix(X), y, alpha = alpha, family = 'binomial', 
                        foldid = folds, intercept = TRUE, type.measure = "auc")

    allobjs[[length(allobjs) + 1]] <- get_best_beta(glmobj)
  }

  bestobj <- get_best_alpha(allobjs, alphalist)

  plot(bestobj$glmobj)
  fit <- bestobj$glmobj$glmnet.fit
  plot(fit, xvar='lambda', label=TRUE)

  return(bestobj)
}

partition_data <- function(outcome, nfolds) {
  require('caret')

  folds <- createFolds(outcome, nfolds)

  foldid <- rep(NA, length(outcome))

  for (ind in 1:nfolds) {
    foldid[folds[[ind]]] <- ind
  }

  return(foldid)
}

get_best_beta <- function(glmobj) {
  # given a glmobj from elastic net, extract the cross-validation score and
  # coefficients of the best fit

  minlambda <- glmobj$lambda.1se
  minlambda.ind <- which(glmobj$lambda == minlambda)
  beta <- fit$beta[, minlambda.ind]
  score <- glmobj$cvm[minlambda.ind]

  return(list(beta = beta, score = score, glmobj = glmobj))
}

get_best_alpha <- function(objlist, alphalist) {
  # given a list of objects returned by get_best_beta and a list of alpha 
  # values, return the object corresponding to the highest score

  scorelist <- sapply(objlist, function(x) {x$score})
  max_ind <- which.max(scorelist)  # if using auc, want max score

  bestobj <- c(objlist[[max_ind]], alpha=alphalist[max_ind])

  return(bestobj)
}