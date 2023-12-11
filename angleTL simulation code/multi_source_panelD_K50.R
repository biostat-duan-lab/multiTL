library(MASS)
library(corpcor)
library(glmnet)

source('~/functions.R')

args <- commandArgs(trailingOnly = TRUE)

M <- as.integer(args[1])
print(paste0('M=', M))

sim = M

simulation <- function(iters, phos_w, vars, p, n){
  mse_beta = mse_w_multi1_lasso = mse_w_multi2 = mse_best_single = 
    mse_multi2 = mse_proposed_list_test = mse_TL_list_test = mse_target_only_list_test = matrix(0, 1, iters)
  
  sim_coef = simulate_coef(var=vars, p=p, rho=phos_w)
  beta = sim_coef$beta
  w = sim_coef$w
  w_unit = apply(w, 2, function(x) x/sqrt(sum(x^2)))
  
  for(j in 1:iters){
    set.seed(sim)
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
    
    mse_beta[,j] = mean_squared_error((x_test %*% beta), y_test)
    
    #best single source
    w_best = w[,which(phos_w==max(phos_w))]
    r = cor(w_best, beta)
    vars = c(var(beta), var(w_best))
    
    fit_best_single = CV_ridge(x_train, x_val, y_train, y_val, w=w_best, rho=r, var=vars)
    beta_best_single = fit_best_single$beta_2D_best
    mse_best_single[,j] = mean_squared_error(x_test%*%beta_best_single, y_test)
    
    #algo 1: obtain w via a validation data
    xy_split = train_test_split(x_train, y_train, size=0.9)
    x_train_1 = data.matrix(xy_split$x_train)
    y_train_1 = data.matrix(xy_split$y_train)
    x_agg = data.matrix(xy_split$x_test)
    y_agg = data.matrix(xy_split$y_test)
    
    y_est = matrix(0, num_K, length(y_agg))
    for(col in 1:num_K){
      y_est[col,] = data.matrix(x_agg) %*% data.matrix(w[,col])
    }

    ### add lasso penalty
    cv.out = cv.glmnet(data.matrix(t(y_est)), y_agg, family = "gaussian", alpha = 0)
    weight_multi1_lasso = cv.out$glmnet.fit$beta[, which(cv.out$lambda == cv.out$lambda.min)]
    weight_multi1_lasso[is.na(weight_multi1_lasso)] = 0
    w_multi1_lasso = w%*%weight_multi1_lasso
    mse_w_multi1_lasso[,j] = mean_squared_error((x_test %*% w_multi1_lasso), y_test)
    
    ##algo 2: ensemble w via eigenvalues
    B = x_train%*%w_unit
    G = cor(B)
    w_weight = (eigen(G)$vectors[,1])^2
    w_multi2 = w_unit%*%w_weight
    mse_w_multi2[,j] = mean_squared_error((x_test %*% w_multi2), y_test)
    
    ##run cross validation using w_multi2
    r = cor(w_multi2, beta)
    vars = c(var(beta), var(w_multi2))
    
    fit_cv_multi2 = CV_ridge(x_train, x_val, y_train, y_val, w=w_multi2, rho=r, var=vars)
    beta_2D_multi2 = fit_cv_multi2$beta_2D_best
    mse_multi2[,j] = mean_squared_error(x_test%*%beta_2D_multi2, y_test)
    
    ##run cross validation using w_multi1_lasso
    r = cor(w_multi1_lasso, beta)
    vars = c(var(beta), var(w_multi1_lasso))
    fit_cv = CV_ridge(x_train_1, x_val, y_train_1, y_val, w=w_multi1_lasso, rho=r, var=vars)
    
    ### angleTL
    beta_2D = fit_cv$beta_2D_best
    mse_proposed_list_test[,j] = mean_squared_error(x_test%*%beta_2D, y_test)
    
    ### distTL:lambda=eta
    beta_1D = fit_cv$beta_1D_best
    mse_TL_list_test[,j] = mean_squared_error(x_test%*%beta_1D, y_test)
    
    ### target-only
    beta_target = fit_cv$beta_target_best
    mse_target_only_list_test[,j] = mean_squared_error(x_test%*%beta_target, y_test)
    
  }
  
  return(list(mse_beta=mse_beta, 
              mse_best_single = mse_best_single,
              mse_w_multi1_lasso=mse_w_multi1_lasso, 
              mse_multi2=mse_multi2,
              mse_proposed_list_test=mse_proposed_list_test, 
              mse_TL_list_test=mse_TL_list_test, 
              mse_target_only_list_test=mse_target_only_list_test))
}


iters = 1
p = 400
n = 2000
vars = c(1,rep(1,5))
num_K = 50
phos_w = c(runif(num_K/2, 0.1, 0.5), 
           runif(num_K/2, 0.6, 0.9))


sim_result = simulation(iters=iters, phos=phos_w, var=vars, p=p, n=n)

save(sim_result, file = paste0('~/results/multi_w_p',p,'_K',num_K,'_sim',sim,".RData"))
