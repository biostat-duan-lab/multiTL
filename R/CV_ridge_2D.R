#' @title CV_ridge_2D
#' @description 3-fold cross validation for selecting tunning parameters
#' @param x_train x for training
#' @param y_train y for training
#' @param w Source model parameter estimators
#' @param rho Correlation of target and source model parametor estimators
#' @param var Variance of target and source model parametor estimators
#' @return A list with the best MSE, lambda and eta
#' @export



CV_ridge_2D <- function(x_train, y_train, w, rho, var){
  mse_2D_best_tmp = list()
  mse_1D_best_tmp = list()
  mse_target_best_tmp = list()

  num=as.integer(nrow(x_train)/3)
  ss=c(rep(1,num),rep(2,num),rep(3,nrow(x_train)-2*num))
  ss=sample(ss)
  n = num*2
  ##############################
  ###Proposed: 2D search
  ##############################
  ##define theoretical optimal lambda & eta
  opt_lam = 1/(n * (1-rho**2))
  opt_eta = rho * opt_lam * sqrt(var[1]/var[2])

  ##define 2D search grid
  lam_2D_grid = c(opt_lam*1000,opt_lam*100, opt_lam*50, opt_lam*20, opt_lam*15, opt_lam*10,opt_lam*5,opt_lam*4,opt_lam*3, opt_lam*2,opt_lam,opt_lam/2,opt_lam/3,opt_lam/4,opt_lam/5, opt_lam/10, opt_lam/15, opt_lam/20, opt_lam/50, opt_lam/100)
  eta_2D_grid = c(opt_eta*1000,opt_eta*100, opt_eta*50, opt_eta*20, opt_eta*15, opt_eta*10,opt_eta*5, opt_eta*4, opt_eta*3, opt_eta*2,opt_eta,opt_eta/2,opt_eta/3,opt_eta/4,opt_eta/5, opt_eta/10, opt_eta/15, opt_eta/20, opt_eta/50, opt_eta/100, 0)

  for (k in 1:3) {

    x = x_train[ss!=k,]
    y = y_train[ss!=k]
    x_val = x_train[ss==k,]
    y_val = y_train[ss==k]

    ##begin 2D grid search
    lam_2D_best = eta_2D_best = beta_2D_best = mse_2D_best = Inf
    mse_2D = matrix(NA,length(lam_2D_grid),length(eta_2D_grid))
    for(i in 1:length(lam_2D_grid)){
      l_2D = lam_2D_grid[i]
      for(j in 1:length(eta_2D_grid)){
        e_2D = eta_2D_grid[j]
        fit_ridge_2D = ridge_closed_form_result(x=x, y=y, lam=l_2D, w=w, eta=e_2D, x_val=x_val, y_val=y_val, standard=FALSE)
        mse_2D[i,j] = fit_ridge_2D$mse
      }
    }
    mse_2D_best_tmp[[k]] = mse_2D
  }

  n = nrow(x_train)
  p = ncol(x_train)

  mse_2D_best0=mse_2D_best_tmp[[1]]+mse_2D_best_tmp[[2]]+mse_2D_best_tmp[[3]]
  lam_2D_best = lam_2D_grid[which(mse_2D_best0==min(mse_2D_best0),arr.ind=TRUE)[1]]
  eta_2D_best = eta_2D_grid[which(mse_2D_best0==min(mse_2D_best0),arr.ind=TRUE)[2]]

  fit_ridge_2D = ridge_closed_form_result(x=x_train, y=y_train, lam=lam_2D_best, w=w, eta=eta_2D_best, x_val=x_train, y_val=y_train, standard=FALSE)
  beta_2D_best = fit_ridge_2D$beta
  mse_2D_best = fit_ridge_2D$mse


  return(list(mse_2D_best=mse_2D_best,beta_2D_best=beta_2D_best,lam_2D_best=lam_2D_best, eta_2D_best=eta_2D_best))
}
