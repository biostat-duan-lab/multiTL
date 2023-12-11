library(MASS)
library(glmnet)
library(pROC)
library(corpcor)

source('~/functions_classification.R')

simulation <- function(iters, phos, vars, p, n){
  auc_beta = auc_w = auc_distTL = auc_angleTL = mse_distTL = mse_angleTL = mis_beta = mis_w = mis_distTL = mis_angleTL = matrix(0, length(phos), iters)
  auc_target = auc_target.glmnet = mis_target.glmnet = mis_target = matrix(0, 1, iters)
  
  for(i in 1:length(phos)){
    set.seed(i)
    r_true = phos[i]
    sim_coef = simulate_coef(vars=vars, p=p, rho=r_true)
    beta = sim_coef$beta
    w_true = sim_coef$w
    for(j in 1:iters){
      #j=1
      print(paste0('rho=',phos[i],'; iteration=',j))
      sim_data = simulate_data(n, p, beta, error=0.1)
      x_train = sim_data$x
      y_train = sim_data$y
      
      sim_val_data = simulate_data(n_val, p, beta, error=0.1)
      x_val = sim_val_data$x
      y_val = sim_val_data$y
      
      sim_test_data = simulate_data(n_test, p, beta, error=0.1)
      x_test = sim_test_data$x
      y_test = sim_test_data$y
      
      sim_src_data = simulate_data(N, p, w_true, error=0.5)
      x_src = sim_src_data$x
      y_src = sim_src_data$y
      
      if(sum(table(y_train)==0)>0 | sum(table(y_val)==0)>0 | sum(table(y_test)==0)>0 | sum(table(y_src)==0)>0){
        next
      }
      
      w = coef(glm(y_src~x_src -1, family = binomial(link='logit')))
      r = cor(w, beta)
      
      ##target only using glmnet
      cv.out = glmnet::cv.glmnet(data.matrix(x_train), y_train, family = "binomial", alpha = 1)
      beta.glmnet = cv.out$glmnet.fit$beta[, which(cv.out$lambda == cv.out$lambda.min)]
      
      ##target only using grid-search
      fit_cv_target = ridge_logistic_1D_CV(x=x_train, y=y_train, beta.init=beta.glmnet, w=NULL, rho=r, vars=NULL, x_val, y_val, standard=TRUE)
      beta_target = ridge_logistic_result(x=x_train, y=y_train, beta.init=beta.glmnet, w=NULL, lambda=fit_cv_target$lam_1D_best, eta=0, standard=TRUE)
      beta_target[is.na(beta_target)]=0
      
      ##distTL
      fit_cv_distTL = ridge_logistic_1D_CV(x=x_train, y=y_train, beta.init=beta_target, w=w, rho=r, vars=vars, x_val, y_val, standard=FALSE)
      beta_distTL = ridge_logistic_result(x=x_train, y=y_train, beta.init=beta_target, w=w, lambda=fit_cv_distTL$lam_1D_best, eta=fit_cv_distTL$lam_1D_best, standard=FALSE)
      beta_distTL[is.na(beta_distTL)]=0
      
      ##angleTL
      fit_cv_angleTL = ridge_logistic_2D_CV(x=x_train, y=y_train, beta.init=beta_target, w=w, rho=r, vars=vars, x_val, y_val)
      beta_angleTL = ridge_logistic_result(x=x_train, y=y_train, beta.init=beta_target, w=w, lambda=fit_cv_angleTL$lam_2D_best, eta=fit_cv_angleTL$eta_2D_best, standard=FALSE)
      beta_angleTL[is.na(beta_angleTL)]=0
      
      ##evaluate using auc as metric
      auc_beta[i,j] = get_auc(invlogit(x_test%*%beta), y_test)
      auc_w[i,j] = get_auc(invlogit(x_test%*%w), y_test)
      auc_target.glmnet[j] = get_auc(invlogit(x_test%*%beta.glmnet), y_test)
      auc_target[j] = get_auc(invlogit(x_test%*%beta_target), y_test)
      auc_distTL[i,j] = get_auc(invlogit(x_test%*%beta_distTL), y_test)
      auc_angleTL[i,j] = get_auc(invlogit(x_test%*%beta_angleTL), y_test)
      
      ##evaluate using misclassification rate
      mis_beta[i,j] = mis_rate(invlogit(x_test%*%beta), y_test)
      mis_w[i,j] = mis_rate(invlogit(x_test%*%w), y_test)
      mis_target.glmnet[j] = mis_rate(invlogit(x_test%*%beta.glmnet), y_test)
      mis_target[j] = mis_rate(invlogit(x_test%*%beta_target), y_test)
      mis_distTL[i,j] = mis_rate(invlogit(x_test%*%beta_distTL), y_test)
      mis_angleTL[i,j] = mis_rate(invlogit(x_test%*%beta_angleTL), y_test)
      
    }
  }
  return(list(auc_beta=auc_beta, auc_w=auc_w, auc_target.glmnet=auc_target.glmnet, auc_target=auc_target, auc_distTL=auc_distTL, auc_angleTL=auc_angleTL,
              mis_beta=mis_beta, mis_w=mis_w, mis_target.glmnet=mis_target.glmnet, mis_target=mis_target, mis_distTL=mis_distTL, mis_angleTL=mis_angleTL))
  
}


iters = 100
phos = c(0.1, 0.3, 0.5, 0.7)
p = 50 #100
n = 50
n_val = 500
n_test = 500
N = 5000
vars = c(1, 15)
vars_name = '1v15'

sim_result = simulation(iters=iters, phos=phos, vars=vars, p=p, n=n)

save(sim_result, file = paste0('~/results/p_',p,'_vars_',vars_name,'_iter_',iters,"_Oct2023.RData"))


