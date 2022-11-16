#' @title mean_squared_error
#' @description Calculating MSE
#' @param y_est,y_test estimated y and actual y
#' @return MSE


mean_squared_error <- function(y_est, y_test){
  mse = mean((y_est - y_test)^2)
  return(mse)
}
