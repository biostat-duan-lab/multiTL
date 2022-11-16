#' @title Agg_fun1
#' @description Aggregating beta estimator using logistic regression
#' @param B /n/home11/yhan/rpackage/multiTL_0.1.0.tar.gz
#' @param x.val,y.val validation X and y
#' @param const A parameter, by default is 2
#' @return Aggregated beta estimator
#' @importFrom stats glm binomial
#' @export

Agg_fun2<-function(B, x.val, y.val, const=2){
  x.val = cbind(1, x.val)
  XX = x.val%*%B
  eta.hat = glm(y.val~XX-1, family = binomial(link = "logit"))$coefficients
  #eta.hat = glm(y.til~XX-1, family = gaussian)$coefficients
  return(eta.hat)
}
