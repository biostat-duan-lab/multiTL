#' @title CV_ridge_target
#' @description 3-fold cross validation for selecting tunning parameters
#' @param x_train x for training
#' @param y_train y for training
#' @return A list with the best MSE, beta and lambda
#' @export



CV_ridge_target <- function(x_train, y_train){
  # mse_2D_best_tmp = list()
  mse_1D_best_tmp = list()
  mse_target_best_tmp = list()
  
  num=as.integer(nrow(x_train)/3)
  ss=c(rep(1,num),rep(2,num),rep(3,nrow(x_train)-2*num))
  ss=sample(ss)
  n = num*2

  ##############################
  ####Traditional ridge, 1D search
  ##############################
  ##define theoretical optimal lambda
  opt_lam_target = 1/n
  
  ##define 1D search grid
  lam_target_grid = c(opt_lam_target*1000,opt_lam_target*100, opt_lam_target*50, opt_lam_target*20, opt_lam_target*15, opt_lam_target*10,opt_lam_target*5,opt_lam_target*4,opt_lam_target*3, opt_lam_target*2,opt_lam_target,opt_lam_target/2,opt_lam_target/3,opt_lam_target/4,opt_lam_target/5, opt_lam_target/10, opt_lam_target/15, opt_lam_target/20, opt_lam_target/50, opt_lam_target/100, opt_lam_target/1000)
  
  for (k in 1:3) {
    
    x = x_train[ss!=k,]
    y = y_train[ss!=k]
    x_val = x_train[ss==k,]
    y_val = y_train[ss==k]
    
    ##begin 1D grid search
    # lam_1D_best = beta_1D_best = mse_1D_best = Inf
    # mse_1D = matrix(NA,length(lam_1D_grid),1)
    # for(i in 1:length(lam_1D_grid)){
    #   l_1D = lam_1D_grid[i]
    #   fit_ridge_1D = ridge_closed_form_result(lam=l_1D, eta=l_1D, x=x, y=y, x_val=x_val, y_val=y_val, w=w, standard=FALSE)
    #   mse_1D[i,1] = fit_ridge_1D$mse
    # }
    # mse_1D_best_tmp[[k]] = mse_1D
    
    ##begin 1D grid search (target only)
    lam_target_best = beta_target_best = mse_target_best = Inf
    mse_target = matrix(NA,length(lam_target_grid),1)
    for(i in 1:length(lam_target_grid)){
      l_target = lam_target_grid[i]
      fit_ridge = ridge_closed_form_result(lam=l_target, eta=NULL, x=x, y=y, x_val=x_val, y_val=y_val, w=NULL, standard=TRUE)
      mse_target[i,1] = fit_ridge$mse
    }
    
    mse_target_best_tmp[[k]] = mse_target
  }
  
  n = nrow(x_train)
  p = ncol(x_train)
  
  # mse_1D_best0=mse_1D_best_tmp[[1]]+mse_1D_best_tmp[[2]]+mse_1D_best_tmp[[3]]
  # lam_1D_best = lam_1D_grid[which(mse_1D_best0==min(mse_1D_best0),arr.ind=TRUE)[1]]
  # 
  # fit_ridge_1D = ridge_closed_form_result(x=x_train, y=y_train, lam=lam_1D_best, w=w, eta=lam_1D_best, x_val=x_train, y_val=y_train, standard=FALSE)
  # beta_1D_best = fit_ridge_1D$beta
  # mse_1D_best = fit_ridge_1D$mse
  
  mse_target_best0=mse_target_best_tmp[[1]]+mse_target_best_tmp[[2]]+mse_target_best_tmp[[3]]
  lam_target_best = lam_target_grid[which(mse_target_best0==min(mse_target_best0),arr.ind=TRUE)[1]]
  
  fit_ridge = ridge_closed_form_result(x=x_train, y=y_train, lam=lam_target_best, w=NULL, eta=NULL, x_val=x_train, y_val=y_train,standard=TRUE)
  beta_target_best = fit_ridge$beta
  mse_target_best = fit_ridge$mse
  
  
  return(list(mse_target_best=mse_target_best, beta_target_best=beta_target_best,lam_target_best=lam_target_best))
}
