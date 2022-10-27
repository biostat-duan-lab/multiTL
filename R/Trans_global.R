#' @title Trans_global
#' @description Derive calibration term estimator and an estimator for the target parameter in each source sites
#' @param X.tar Variables from the target. The variables need to be completely the same set and in the same order as variables used in the model parameter estimators.
#' @param y.tar Response from the target.
#' @param X.src A list of shareable source X, the length of the list should be M
#' @param y.src A list of shareable source y, the length of the list should be M
#' @param delta A vector of calibration estimator for each source
#' @importFrom glmnet cv.glmnet
#' @return A list of source estimator, calibration term estimator and an estimator for the target parameter
#' @export

Trans_global<-function(X.tar, y.tar, X.src, y.src, delta=NULL){
  p <- ncol(X.tar)
  n0.tar <- length(y.tar)
  K <- length(y.src)  ###(K>1)
  ##global method
  Xdelta = c()
  if(is.null(delta)){
    for(k in 1:K){
      w.k = as.numeric(ST_init(X.src[[k]], y.src[[k]])$beta0)
      fit.tar <- cv.glmnet(x=X.tar, y=y.tar, nfolds=5, family='binomial', offset=w.k[1]+X.tar%*%w.k[-1], lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar))
      delta.k = c(fit.tar$glmnet.fit$a0[which(fit.tar$lambda == fit.tar$lambda.min)], fit.tar$glmnet.fit$beta[, which(fit.tar$lambda == fit.tar$lambda.min)])
      delta.k.thre = thres(delta.k, sqrt(n0.tar), p) ###+threshold
      Xdelta.k = tcrossprod(delta.k.thre[-1], X.src[[k]])+delta.k.thre[1]
      Xdelta = c(Xdelta, Xdelta.k)
    }
  }else{
    for(k in 1:K){
      Xdelta.k = tcrossprod(delta[[k]][-1], X.src[[k]])+delta[[k]][1]
      Xdelta = c(Xdelta, Xdelta.k)
    }
  }

  XX.src <- yy.src <- NULL
  for(k in 1:K){
    XX.src <- rbind(XX.src, X.src[[k]])
    yy.src <- c(yy.src, y.src[[k]])
  }

  n0.tar.global <- length(c(y.tar,yy.src))
  offset <- c(rep(0, nrow(X.tar)), -Xdelta)
  fit.global <- cv.glmnet(x=rbind(X.tar,XX.src), y=c(y.tar,yy.src), nfolds=5, family='binomial', offset=offset, lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar.global))
  beta.hat = c(fit.global$glmnet.fit$a0[which(fit.global$lambda == fit.global$lambda.min)], fit.global$glmnet.fit$beta[, which(fit.global$lambda == fit.global$lambda.min)])

  return(beta.hat)
}
