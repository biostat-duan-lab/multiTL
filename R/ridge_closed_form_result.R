#' @title ridge_closed_form_result
#' @description Calculating closed form result in a ridge regression
#' @param x x for training
#' @param y y for training
#' @param x_val x for validation
#' @param y_val y for validation
#' @param lam,eta Tuning parameter
#' @param w Source model parameter estimates
#' @param standard TRUE: run traditional ridge regression
#' @return A list with the beta and MSE
#' @export


ridge_closed_form_result <- function(x, y, lam, w, eta, x_val, y_val, standard){
  n = nrow(x)
  p = ncol(x)
  I = diag(p)


  if(n>p){
    # print(dim(t(x)%*%x + n * lam * I))

    if(standard==TRUE){
      beta_est = solve(t(x)%*%x + n * lam * I, tol=1e-1000) %*% t(x) %*% y
    }else{
      w_vector = matrix(w, p, 1)
      beta_est = solve(t(x)%*%x + n * lam * I, tol=1e-1000) %*% (t(x) %*% y + n * eta * w_vector)
    }
  }else{

    x_svd = svd(x, nu=n, nv=p)
    u = x_svd$u
    s = cbind(diag(x_svd$d), matrix(0, n, p-n))
    v = x_svd$v
    R = u %*% s
    RtR = t(R) %*% R
    RtY = t(R) %*% y

    if(corpcor::is.positive.definite(RtR + n * lam * I)==FALSE){
      cov = corpcor::make.positive.definite(RtR + n * lam * I)
    }else{
      cov = RtR + n * lam * I
    }
    if(standard==TRUE){
      beta_est = v %*% solve(cov, tol=1e-1000) %*% RtY
    }else{
      w_vector = matrix(w, p, 1)
      beta_est = v %*% solve(cov, tol=1e-1000) %*% (RtY + n * eta * t(v) %*% w_vector)
    }

  }
  mse = mean_squared_error(x_val%*%beta_est, y_val)
  return(list(beta=beta_est, mse=mse))
}
