#' @title target_only
#' @description target_only
#' @param X Variables from the target. The variables need to be completely the same set and in the same order as variables used in the model parameter estimators.
#' @param y Response from the target.
#' @return A list of effect estimator and tuning parameters lambda and eta from target-only.
#' @importFrom stats coef cor lm
#' @export

target_only <- function(X, y){
  
  glm.tar = glmnet(X,y,family = 'gaussian',alpha=0)
  beta <- predict(glm.tar, s=0.05, type = 'coefficients')[-1]
  
  mse_w = mse_proposed_list = mse_proposed_list_test = mse_TL_list = mse_TL_list_test = NA
  lam_proposed_list = eta_proposed_list = lam_TL_list = NA
  mse_target_only_list = mse_target_only_list_test = lam_target_list = NA
  
  xy_split = train_test_split(X, y, size=0.6)
  x_train = xy_split$x_train
  y_train = xy_split$y_train
  x_val = xy_split$x_test
  y_val = xy_split$y_test
  
  ##run cross validation
  fit_cv = CV_ridge_target(x_train=x_train, y_train=y_train)
  mse_target_best = fit_cv$mse_target_best
  beta_target_best = fit_cv$beta_target_best
  lam_target_best = fit_cv$lam_target_best
  
  return(list(beta=beta_target_best,lam=lam_target_best))
  
}
