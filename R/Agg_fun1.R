#' @title Agg_fun1
#' @description Aggregating beta estimator
#' @param B A matrix of target-only estimator and target estimator from each source in each column
#' @param x.val,y.val validation X and y
#' @param const A parameter, by default is 2
#' @return Aggregated beta estimator
#' @export

Agg_fun1<-function(B, x.val, y.val, const=2){
  x.val = cbind(1, x.val)
  loss.B <- apply(B, 2, function(b) - sum(y.val*log(logistic(x.val%*%b))+(1-y.val)*log(1-logistic(x.val%*%b))))
  eta.hat <- exp(-const*loss.B)/sum(exp(-const*loss.B))
  return(eta.hat)
}
