#' @title angleTL
#' @description Angle-based transfer learning
#' @param X Variables from the target. The variables need to be completely the same set and in the same order as variables used in the model parameter estimators.
#' @param y Response from the target.
#' @param w.src Pre-trained source model parameter estimators.
#' @return A list of effect estimator and tuning parameters lambda and eta from angleTL.
#' @importFrom stats coef cor lm
#' @importFrom glmnet glmnet
#' @export

angleTL <- function(X, y, w.src){
  glm.tar = glmnet(X,y,family = 'gaussian',alpha=0)
  beta <- predict(glm.tar, s=0.05, type = 'coefficients')[-1]
  var = c(var(beta),var(w.src))
  rho = cor(beta, w.src)
  mse_w = mse_proposed_list = mse_proposed_list_test = mse_TL_list = mse_TL_list_test = NA
  lam_proposed_list = eta_proposed_list = lam_TL_list = NA
  mse_target_only_list = mse_target_only_list_test = lam_target_list = NA

  xy_split = train_test_split(X, y, size=0.6)
  x_train = xy_split$x_train
  y_train = xy_split$y_train
  x_val = xy_split$x_test
  y_val = xy_split$y_test

  ##run cross validation
  fit_cv = CV_ridge_2D(x_train=x_train, y_train=y_train,w=w.src, rho=rho, var=var)
  mse_2D_best = fit_cv$mse_2D_best
  beta_2D_best = fit_cv$beta_2D_best
  lam_2D_best = fit_cv$lam_2D_best
  eta_2D_best = fit_cv$eta_2D_best

  return(list(beta=beta_2D_best,lam=lam_2D_best, eta=eta_2D_best ))

}
