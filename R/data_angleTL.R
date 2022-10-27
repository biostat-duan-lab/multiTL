#' Simulated data for angleTL
#' @name data_angleTL
#' @docType data
#' @format A list with target x and y, source estimator
NULL

# library(MASS)
# library(corpcor)
# library(glmnet)
#
# simulate_coef <- function(mean=c(0,0), var, rho, p,ortho = F){
#   if(ortho){
#     sd_delta = sqrt((1/(rho^2)-1)*var[2])
#     beta = rnorm(p,mean = 0,sd = sqrt(var[1]))
#     # w = rnorm(p,mean = 0,sd = sqrt(var[2]))
#     delta = rnorm(p, mean=0, sd = sd_delta)
#     w = beta - delta
#   }
#   else{
#     cov_raw = matrix(rbind(c(var[1], sqrt(var[1]*var[2])*rho),c(sqrt(var[1]*var[2])*rho, var[2])), 2, 2)
#     if(corpcor::is.positive.definite(cov_raw)==FALSE){
#       cov = corpcor::make.positive.definite(cov_raw)
#     }else{
#       cov = cov_raw
#     }
#     beta_w = matrix(mvrnorm(n=p, mu=mean, Sigma=cov), p, 2)
#     w = matrix(beta_w[,2],p,1)
#     beta =  matrix(beta_w[,1],p,1)
#   }
#   return(return(list(beta=beta, w=w)))
# }
#
# var=c(1,1)
# rho=0.6
# p=50
#
# coef = simulate_coef(mean=c(0,0), var, rho, p,ortho = F)
# w.src=coef$w
# beta.tar=coef$beta
#
# simulate_data <- function(n, p, beta, error=0.5){
#   x = matrix(rnorm(n*p, 0, 1), n, p)
#   y = matrix(rnorm(n, x%*%beta, error), n, 1)
#   return(list(x=x, y=y))
# }
#
# n=50
# data = simulate_data(n, p, beta.tar, error=0.5)
# X=data$x
# y=data$y
# w.src=coef$w
# data_angleTL=list(X=X,y=y,w.src=w.src)

