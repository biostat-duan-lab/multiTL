#' @title ST_init
#' @description Derive beta estimator from target data
#' @param y.tar Target y
#' @param X.tar Target X
#' @importFrom glmnet cv.glmnet
#' @return A list of beta and tuning parameter lambda
#' @export

ST_init<-function(X.tar,y.tar){
  #single-task initialization
  p <- ncol(X.tar)
  n0.tar <- length(y.tar)
  fit <- cv.glmnet(x=X.tar, y=y.tar, nfolds = 4, family='binomial')
  lam.const = fit$lambda.min / sqrt(2*log(p)/n0.tar)
  beta0 = c(fit$glmnet.fit$a0[which(fit$lambda == fit$lambda.min)], fit$glmnet.fit$beta[, which(fit$lambda == fit$lambda.min)])
  return(list(beta0=as.numeric(beta0), lam.const=as.numeric(lam.const)))
}
