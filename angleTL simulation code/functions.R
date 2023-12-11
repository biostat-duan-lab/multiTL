mean_squared_error <- function(y_est, y_test){
  mse = mean((y_est - y_test)^2)
  return(mse)
}

train_test_split <- function(x, y, size){
  n = length(y)
  train = sample(1:n, size*n)
  test = 1:n
  test = test[!test %in% train]
  x_train = x[train,]
  y_train = y[train]
  x_test = x[test,]
  y_test = y[test]
  return(list(x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test))
}

# simulate_coef <- function(mean=0, var, rho, p){
#   beta = matrix(mvrnorm(n=p, mu=mean, Sigma=var[1]), p, 1)
#   sd_delta = sqrt(1/(rho^2)-1)
#   
#   delta = rnorm(p, mean=mean, sd = sd_delta)
#   w = beta + delta
#   
#   return(list(beta=beta, w=w))
# }

simulate_coef <- function(mean=c(0,0), vars, rho, p){
  cov_raw = matrix(rbind(c(vars[1], sqrt(vars[1]*vars[2])*rho),c(sqrt(vars[1]*vars[2])*rho, vars[2])), 2, 2)
  if(corpcor::is.positive.definite(cov_raw)==FALSE){
    cov = corpcor::make.positive.definite(cov_raw)
  }else{
    cov = cov_raw
  }
  beta_w = matrix(mvrnorm(n=p, mu=mean, Sigma=cov), p, 2)
  
  return(list(beta=matrix(beta_w[,1],p,1), w=matrix(beta_w[,2],p,1)))
}


ar1_cor <- function(n, rho){
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}

simulate_data <- function(n, p, beta, error=0.5){
  sigma = ar1_cor(p, 0.3)
  x = matrix(rnorm(n*p, 0, sigma), n, p)
  y = matrix(rnorm(n, x%*%beta, error), n, 1)
  y_center = (y - mean(y))/sd(y)
  return(list(x=x, y=y_center))
}


ridge_closed_form_result <- function(n, p, v, RtR, lam, I, RtY, w_vector, eta, x_val, y_val, standard){
  if(standard==TRUE){
    beta_est = v %*% solve(RtR + n * lam * I, tol=1e-100) %*% RtY
  }else{
    beta_est = v %*% solve(RtR + n * lam * I, tol=1e-100) %*% (RtY + n * eta * t(v) %*% w_vector)
  }
  mse = mean_squared_error(x_val%*%beta_est, y_val)
  
  return(list(beta=beta_est, mse=mse))
}

CV_ridge <- function(x_train, x_val, y_train, y_val, w, rho, var){
  
  ##############################
  ###Proposed: 2D search
  ##############################
  n = nrow(x_train)
  p = ncol(x_train)
  
  x_svd = svd(x_train, nu=n, nv=p)
  u = x_svd$u
  s = cbind(diag(x_svd$d), matrix(0, n, p-n))
  v = x_svd$v
  w_vector = matrix(w, p, 1)
  R = u %*% s
  RtR = t(R) %*% R
  RtY = t(R) %*% y_train
  I = diag(p)
  
  ##define theoretical optimal lambda & eta
  opt_lam = 1/(n * var[1] *(1-rho**2))
  #opt_eta = rho * opt_lam * sqrt(var[1]/var[2])
  
  ##define 2D search grid
  lam_2D_grid = c(opt_lam*10, opt_lam*5, opt_lam*2, opt_lam*1.5, opt_lam, opt_lam/1.5, opt_lam/2, opt_lam/5, opt_lam/10, opt_lam/100, seq(0.0001, 0.5, 0.005))
  #eta_2D_grid = c(opt_eta*10, opt_eta*5, opt_eta*2, opt_eta*1.5, opt_eta, opt_eta/1.5, opt_eta/2, opt_eta/5, opt_eta/10, opt_eta/100)
  
  ##begin 2D grid search
  lam_2D_best = eta_2D_best = beta_2D_best = mse_2D_best = Inf
  for(i in 1:length(lam_2D_grid)){
    l_2D = lam_2D_grid[i]
    eta_2D_grid = c(rho*l_2D*sqrt(var[1]/var[2])*5, rho*l_2D*sqrt(var[1]/var[2])*3, rho*l_2D*sqrt(var[1]/var[2])*2, rho*l_2D*sqrt(var[1]/var[2]), rho*l_2D*sqrt(var[1]/var[2])/2, rho*l_2D*sqrt(var[1]/var[2])/3, rho*l_2D*sqrt(var[1]/var[2])/5, 0)
    for(j in 1:length(eta_2D_grid)){
      e_2D = eta_2D_grid[j]
      fit_ridge_2D = ridge_closed_form_result(n, p, v, RtR, lam=l_2D, I, RtY, w_vector, eta=e_2D, x_val=x_val, y_val=y_val, standard=FALSE)
      beta_2D = fit_ridge_2D$beta
      mse_2D = fit_ridge_2D$mse
      if(mse_2D < mse_2D_best){
        lam_2D_best = l_2D
        eta_2D_best = e_2D
        beta_2D_best = beta_2D
        mse_2D_best = mse_2D 
      }
    }
  }
  
  ##############################
  ####TL: lambda=eta, 1D search 
  ##############################
  ##define theoretical optimal lambda
  opt_lam_1D = 1/(n * var[1] *(1-rho**2))
  
  ##define 1D search grid
  lam_1D_grid = c(opt_lam_1D, opt_lam_1D*2, opt_lam_1D*5, opt_lam_1D*10, opt_lam_1D*100, opt_lam_1D/2, opt_lam_1D/5, opt_lam_1D/10, opt_lam_1D/100, opt_lam_1D/10000, seq(0.0001, 0.5, 0.005))
  
  ##begin 1D grid search
  lam_1D_best = beta_1D_best = mse_1D_best = Inf
  for(i in 1:length(lam_1D_grid)){
    l_1D = lam_1D_grid[i]
    fit_ridge_1D = ridge_closed_form_result(n, p, v, RtR, lam=l_1D, I, RtY, w_vector, eta=l_1D, x_val=x_val, y_val=y_val, standard=FALSE)
    beta_1D = fit_ridge_1D$beta
    mse_1D = fit_ridge_1D$mse
    if(mse_1D < mse_1D_best){
      lam_1D_best = l_1D
      mse_1D_best = mse_1D
      beta_1D_best = beta_1D
    }
  }
  
  ##############################
  ####Traditional ridge
  ##############################
  ##define theoretical optimal lambda
  opt_lam_target = 1/n
  
  ##define 1D search grid
  lam_target_grid = c(seq(0.0001, 0.5, 0.005), opt_lam_target*10,opt_lam_target*5, opt_lam_target*2, opt_lam_target, opt_lam_target/100,opt_lam_target/10,opt_lam_target/2,opt_lam_target/5)
  
  ##begin 1D grid search
  lam_target_best = beta_target_best = mse_target_best = Inf
  for(i in 1:length(lam_target_grid)){
    l_target = lam_target_grid[i]
    fit_ridge = ridge_closed_form_result(n, p, v, RtR, lam=l_target, I, RtY, w_vector=NULL, eta=NULL, x_val=x_val, y_val=y_val, standard=TRUE)
    beta_target = fit_ridge$beta
    mse_target = fit_ridge$mse
    if(mse_target < mse_target_best){
      lam_target_best = l_target
      beta_target_best = beta_target
      mse_target_best = mse_target
    }
  }
  
  return(list(mse_2D_best=mse_2D_best, mse_1D_best=mse_1D_best, mse_target_best=mse_target_best, beta_2D_best=beta_2D_best, beta_1D_best=beta_1D_best, beta_target_best=beta_target_best))
}
