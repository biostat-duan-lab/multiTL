#' @title logistic
#' @description Get the logistic form of x
#' @param x Original x
#' @return logistic x

logistic<-function(x){
  1/(1+exp(-x))
}
