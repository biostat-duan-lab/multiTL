#' @title TL_init
#' @description Derive calibration term estimator and an estimator for the target parameter in each source sites
#' @param y.tar Target y
#' @param X.tar Target X
#' @param y.src Source y if available, by default is NULL
#' @param X.src Source X if available, by default is NULL
#' @param w Source estimator if provided, by default is NULL
#' @importFrom glmnet cv.glmnet
#' @return A list of source estimator, calibration term estimator and an estimator for the target parameter
#' @export
#'
TL_init<-function(X.tar, y.tar, X.src=NULL, y.src=NULL, w=NULL){
  p <- ncol(X.tar)
  n0.src <- length(y.src)
  n0.tar <- length(y.tar)

  if(is.null(w)){
    #source population estimates
    fit.src <- cv.glmnet(x=X.src, y=y.src, family='binomial', nfolds=4, lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar))
    lam.const = fit.src$lambda.min / sqrt(2*log(p)/n0.src)
    w0 <- c(fit.src$glmnet.fit$a0[which(fit.src$lambda == fit.src$lambda.min)], fit.src$glmnet.fit$beta[, which(fit.src$lambda == fit.src$lambda.min)])
  }else{
    w0 = w
  }

  #target population estimates
  fit.tar <- cv.glmnet(x=X.tar, y=y.tar, nfolds=5, family='binomial', offset=w0[1]+X.tar%*%w0[-1], lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar))
  delta0 = c(fit.tar$glmnet.fit$a0[which(fit.tar$lambda == fit.tar$lambda.min)], fit.tar$glmnet.fit$beta[, which(fit.tar$lambda == fit.tar$lambda.min)])

  return(list(w0=w0, delta0=delta0, beta0=w0+delta0))
}
