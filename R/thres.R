#' @title thres
#' @description Shrink estimators of somr sites with performance lower than a certain threshold to zero
#' @param b Original calibration estimator
#' @param k Target sample size
#' @param p Dimension
#' @return Calibration estimator after adding threshold
#' @export

thres<-function(b, k, p){
  b*(rank(abs(b))>=(p-k))
}
