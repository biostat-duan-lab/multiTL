
mean_squared_error <- function(y_est, y_test){
  mse = mean((y_est - y_test)^2)
  return(mse)
}

get_r2 <- function(y_est, y_test){
  r2 = (cor(y_est, y_test))^2
  return(r2)
}

get_auc <- function(y_hat, y_test){
  auc = auc(pROC::roc(y_test, as.numeric(y_hat), quiet=TRUE))
  return(as.numeric(auc))
}

mis_rate <- function(p_hat, y_test){
  y_hat = ifelse(p_hat > 0.5, 1, 0)
  misclassification_rate = sum(y_hat != y_test) / length(y_test)
  return(as.numeric(misclassification_rate))
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


invlogit <- function(x){
  1 / (1 + exp(-(x)))
}

simulate_data <- function(n, p, beta, error=0.5){
  x = matrix(rnorm(n*p, 0, 1), n, p)
  prob_y = invlogit(x%*%beta)
  y = matrix(rbinom(n, 1, prob_y), n, 1)
  x_center = apply(x, 2, function(x) x/sqrt(sum(x^2)))
  y_center = y
  return(list(x=x_center, y=y_center))
}



ridge_logistic_result <- function(x, y, beta.init, w, lambda, eta, standard=FALSE){
  beta = beta.init
  if(is.null(w)==FALSE){w_vector = matrix(w, dim(x)[2], 1)}
  
  # Begin Newton-Raphson loop
  for (rep in 1:5) {
    
    # Calculate predicted probabilities
    pr <- exp(x %*% beta) / (1 + exp(x %*% beta))
    
    # Calculate the diagonal weight matrix
    W <- diag(as.vector(pr * (1 - pr)))
    
    # Calculate the Hessian & gradient
    if(standard==TRUE){
      H <- t(x) %*% W %*% x
      grad <- t(x) %*% (y - pr)
    }else{
      H <- t(x) %*% W %*% x + 2 * lambda * diag(ncol(x))
      grad <- t(x) %*% (y - pr) - 2 * lambda * beta + 2 * eta * w_vector
    }
    
    # Update beta
    diff = solve(H, grad, tol=1e-100)
    beta <- beta + diff
    #print(mean(diff))
    rep = rep + 1
  }
  
  beta_est = beta
  
  return(beta_est)
} 


ridge_logistic_2D_CV <- function(x, y, beta.init, w, rho, vars, x_val, y_val){
  n = length(y)
  opt_lam = 1/(n * vars[1] * (1-rho**2))
  opt_eta = rho * opt_lam * sqrt(vars[1]/vars[2])
  
  lam_2D_grid = c(seq(0.0001, 0.5, 0.005), 
                  opt_lam*100, opt_lam*80, opt_lam*50, opt_lam*20, opt_lam*10, 
                  seq(0,5,0.1)*opt_lam,
                  opt_lam/10, opt_lam/20, opt_lam/50, opt_lam/80, opt_lam/100, opt_lam/1000)
  eta_2D_grid = c(rho * lam_2D_grid * sqrt(vars[1]/vars[2]),
                  lam_2D_grid, 0)
  
  ##begin 2D grid search
  ##find the optimal lambda, given 2 eta values
  lam_2D_best = eta_2D_best = beta_2D_best = auc_2D_best = -Inf
  for(i in 1:length(lam_2D_grid)){
    l_2D = lam_2D_grid[i]
    for(j in 1:length(eta_2D_grid)){
      e_2D = eta_2D_grid[j]
      beta_2D = ridge_logistic_result(x, y, beta.init, w, lambda=l_2D, eta=e_2D, standard=FALSE)
      beta_2D[is.na(beta_2D)] = 0
      auc_2D = get_auc(invlogit(x_val%*%beta_2D), y_val)
      #print(paste0('current auc=',auc_2D))
      #print(paste0('best auc=',auc_2D_best))
      #print(paste0('opt lam=', lam_2D_best,'; opt eta=',eta_2D_best))
      if(auc_2D > auc_2D_best){
        lam_2D_best = l_2D
        eta_2D_best = e_2D
        beta_2D_best = beta_2D
        auc_2D_best = auc_2D 
      }
    }
  }
  
  # ##search for the optimal eta, given the optimal lambda
  # opt_eta = rho*lam_2D_best*sqrt(vars[1]/vars[2])
  # eta_extra = c(opt_eta*100, opt_eta*10, opt_eta*5, opt_eta*2, opt_eta/2, opt_eta/5, opt_eta/10, opt_eta/100, lam_2D_best, lam_2D_best/2, lam_2D_best*2, lam_2D_best*10, lam_2D_best/10)
  # for(k in 1:length(eta_extra)){
  #     e_2D = eta_extra[k]
  #     beta_2D = ridge_logistic_result(x, y, beta.init, w, lambda=lam_2D_best, eta=e_2D, standard=FALSE)
  #     beta_2D[is.na(beta_2D)] = 0
  #     auc_2D = get_auc(invlogit(x_val%*%beta_2D), y_val)
  #     print(paste0('current auc=',auc_2D))
  #     print(paste0('best auc=',auc_2D_best))
  #     print(paste0('opt lam=', lam_2D_best,'; opt eta=',eta_2D_best))
  #     if(auc_2D > auc_2D_best){
  #         eta_2D_best = e_2D
  #         beta_2D_best = beta_2D
  #         auc_2D_best = auc_2D
  #     }
  # }
  
  return(list(auc_2D_best=auc_2D_best, beta_2D_best=beta_2D_best, lam_2D_best=lam_2D_best, eta_2D_best=eta_2D_best))
}




ridge_logistic_1D_CV <- function(x, y, beta.init, w, rho, vars, x_val, y_val, standard){
  n = length(y)
  if(standard==FALSE){
    opt_lam = 1/(n * vars[1] * (1-rho**2))
  }else{
    opt_lam = 1/n
  }
  lam_1D_grid = c(seq(0.0001, 0.5, 0.005), 
                  opt_lam*100, opt_lam*80, opt_lam*50, opt_lam*20, opt_lam*10, 
                  seq(0,5,0.1)*opt_lam,
                  opt_lam/10, opt_lam/20, opt_lam/50, opt_lam/80, opt_lam/100, opt_lam/1000)
  
  ##begin 1D grid search
  ##find the optimal lambda that maximize the auc in a validation data
  lam_1D_best = beta_1D_best = auc_1D_best = -Inf
  for(i in 1:length(lam_1D_grid)){
    l_1D = lam_1D_grid[i]
    if(standard==FALSE){
      e_1D = l_1D
    }else{
      e_1D = 0
    }
    beta_1D = ridge_logistic_result(x, y, beta.init, w, lambda=l_1D, eta=e_1D, standard)
    beta_1D[is.na(beta_1D)] = 0
    auc_1D = get_auc(invlogit(x_val%*%beta_1D), y_val)
    #print(paste0('current auc=',auc_1D))
    #print(paste0('best auc=',auc_1D_best))
    #print(paste0('opt lam=', lam_1D_best))
    if(auc_1D > auc_1D_best){
      lam_1D_best = l_1D
      beta_1D_best = beta_1D
      auc_1D_best = auc_1D 
    }
  }
  return(list(auc_1D_best=auc_1D_best, beta_1D_best=beta_1D_best, lam_1D_best=lam_1D_best))
}