
dolfpglm <- function(dfile){
  # given an input data file, perform an elastic net regression
  # on the data contained therein
  dat <- read.table(dfile, sep=',', header=TRUE, row.names=1)
  
  # make model matrix
  y <- dat[, 1]
  X <- dat[, -1]
  
  # run sparse regression
  # don't include constant in X; glmnet will fit an intercept, but if we supply
  # one manually, it will be penalized, which we don't want
  glmobj <- cv.glmnet(as.matrix(X), y, family = 'binomial', 
                      nfolds=10, intercept=TRUE, type.measure="auc")
  plot(glmobj)
  fit <- glmobj$glmnet.fit
  plot(fit, xvar='lambda', label=TRUE)
  minlambda <- glmobj$lambda.min
  minlambda.ind <- which(glmobj$lambda == minlambda)
  beta <- fit$beta[, minlambda.ind]
  return (list(beta = beta, glmobj = glmobj))
}
