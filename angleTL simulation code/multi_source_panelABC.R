library(MASS)
library(corpcor)
library(glmnet)

source('~/functions.R')

args <- commandArgs(trailingOnly = TRUE)

M <- as.integer(args[1])
print(paste0('M=', M))

sim = M

simulation <- function(iters, phos_w, vars, p, n, N){
  mse_beta = mse_w1 = mse_w2 = mse_w3 = mse_w4 = mse_w5 = mse_w_naive_ensemble = mse_w_ensemble = 
    mse_best_single = mse_multi1 = mse_multi2 = mse_distTL = mse_target = matrix(0, 1, iters)
  
  sim_coef = simulate_coef(var=vars, p=p, rho=phos_w)
  beta = sim_coef$beta
  w = sim_coef$w

  for(j in 1:iters){
    print(j)
    ##generate w estimates
    set.seed(j)
    w_est = matrix(0,p,5)
    for(k in 1:5){
      sim_data = simulate_data(N, p, w[,k])
      x_w = sim_data$x
      y_w = sim_data$y
      lm_w = lm(y_w~x_w) 
      w_est[,k] = coef(lm_w)[-1]
    }
    w_est[is.na(w_est)]=0
    w_unit = apply(w_est, 2, function(x) x/sqrt(sum(x^2)))
    
    ##generate target data
    sim_data = simulate_data(n, p, beta)
    x_norm = sim_data$x
    y_norm = sim_data$y
    
    xy_split = train_test_split(x_norm, y_norm, size=0.6)
    x_train_tmp = xy_split$x_train
    y_train_tmp = xy_split$y_train
    x_val = xy_split$x_test
    y_val = xy_split$y_test
    
    xytrain_split = train_test_split(x_train_tmp, y_train_tmp, size=0.2)
    x_train = xytrain_split$x_train
    y_train = xytrain_split$y_train
    x_test = xytrain_split$x_test
    y_test = xytrain_split$y_test
    
    ##MSE of source-only
    mse_beta[,j] = mean_squared_error((x_test %*% beta), y_test)
    mse_w1[,j] = mean_squared_error((x_test %*% w_unit[,1]), y_test)
    mse_w2[,j] = mean_squared_error((x_test %*% w_unit[,2]), y_test)
    mse_w3[,j] = mean_squared_error((x_test %*% w_unit[,3]), y_test)
    mse_w4[,j] = mean_squared_error((x_test %*% w_unit[,4]), y_test)
    mse_w5[,j] = mean_squared_error((x_test %*% w_unit[,5]), y_test)
    
    ##algo 1: naive ensemble w via training data
    loss.w = apply(w, 2, function(w) 1/sum((y_train-x_train%*%w)^2)) #####correlation(y_hat, y_train) to select weights, increase # of w, and correlation
    w_naive_weight = loss.w/sum(loss.w)
    w_naive_ensemble = w_unit%*%w_naive_weight
    mse_w_naive_ensemble[,j] = mean_squared_error((x_test %*% w_naive_ensemble), y_test)
    
    ##algo 2: ensemble w via eigenvalues
    B = x_train%*%w_unit
    G = cor(B)
    w_weight = (eigen(G)$vectors[,1])^2
    w_ensemble = w_unit%*%w_weight
    mse_w_ensemble[,j] = mean_squared_error((x_test %*% w_ensemble), y_test)
    
    r = cor(w_ensemble, beta)
    vars = c(var(beta), var(w_ensemble))
    
    ##run cross validation using the best single source w
    fit_cv = CV_ridge(x_train, x_val, y_train, y_val, w=w_unit[,5], rho=r, var=vars)
    beta_2D_best = fit_cv$beta_2D_best
    mse_best_single[,j] = mean_squared_error(x_test%*%beta_2D_best, y_test)
    
    ##run cross validation using naive ensemble w via training data
    fit_cv = CV_ridge(x_train, x_val, y_train, y_val, w=w_naive_ensemble, rho=r, var=vars)
    beta_2D_best = fit_cv$beta_2D_best
    mse_multi1[,j] = mean_squared_error(x_test%*%beta_2D_best, y_test)
    
    ##run cross validation using w_ensemble
    fit_cv = CV_ridge(x_train, x_val, y_train, y_val, w=w_ensemble, rho=r, var=vars)
    
    ### angleTL
    beta_2D = fit_cv$beta_2D_best
    mse_multi2[,j] = mean_squared_error(x_test%*%beta_2D, y_test)
    
    ### distTL:lambda=eta
    beta_1D = fit_cv$beta_1D_best
    mse_distTL[,j] = mean_squared_error(x_test%*%beta_1D, y_test)
    
    ### target-only
    beta_target = fit_cv$beta_target_best
    mse_target[,j] = mean_squared_error(x_test%*%beta_target, y_test)
    
  }
  return(list(mse_beta=mse_beta, 
              mse_w1=mse_w1, mse_w2=mse_w2, mse_w3=mse_w3, mse_w4=mse_w4, mse_w5=mse_w5, 
              mse_w_naive_ensemble=mse_w_naive_ensemble, mse_w_ensemble=mse_w_ensemble,
              mse_best_single=mse_best_single,
              mse_multi1=mse_multi1, 
              mse_multi2=mse_multi2, 
              mse_distTL=mse_distTL, 
              mse_target=mse_target))
}



iters = 100
p = 50
n = 100
N = 5000
vars = c(1,rep(1,5))

####panel A
phos_w = c(0.4, 0.45, 0.5, 0.55, 0.6)
r_name = '0.4v0.45v0.5v0.55v0.6'

####panel B
phos_w = c(0.1, 0.3, 0.5, 0.7, 0.9)
r_name = '0.1v0.3v0.5v0.7v0.9'

####panel C
phos_w = c(0.1,0.1,0.1,0.1,0.1)
r_name = '0.1v0.1v0.1v0.1v0.1'

sim_result = simulation(iters=iters, phos=phos_w, var=vars, p=p, n=n, N=N)

save(sim_result, file = paste0('~/multi_w_p_',p,'_r_',r_name,'_iter_',iters,'_N_',N,"_noise.RData"))


