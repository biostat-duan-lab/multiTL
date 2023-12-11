library(MASS)
library(pracma)
library(corpcor)
library(grDevices)
library(ggplot2)
library(patchwork)

import_path = '/Users/tiangu/OneDrive - Harvard University/Tian Gu-Shared/Angle-based Ridge/'

simulate_ridge_risk_cov = function(Sigma, n, p, rho, lambda_arr, sd_delta, Vt=1, Vs=1, sigma2=1){
  
  ## Monte Carlo Evaluation
  num_monte = 500      #5*1e2
  
  a2_t = p*Vt
  a2_s = p*Vs
  gamma = p/n
  S = Sigma^(1/2)
  
  pred_risk_t = pred_risk_1D = pred_risk_2D = zeros(length(lambda_arr),1)
  for(k in 1:length(lambda_arr)){
    print(k)
    lambda = lambda_arr[k]
    pred_err_t = pred_err_1D = pred_err_2D = zeros(num_monte,1)
    for(i in 1:num_monte){
      X = pracma::randn(n,p)%*%S
      # beta = sqrt(Vt) * pracma::randn(p,1)
      # #generate w
      # sd_delta = sqrt(Vt/(rho^2)-Vt)
      # delta = rnorm(p, mean=0, sd = sd_delta)
      # w = beta + delta
      
      cov_raw = matrix(rbind(c(Vt, sqrt(Vt*Vs)*rho),c(sqrt(Vt*Vs)*rho, Vs)), 2, 2)
      if(corpcor::is.positive.definite(cov_raw)==FALSE){
        cov = corpcor::make.positive.definite(cov_raw)
      }else{
        cov = cov_raw
      }
      beta_w = matrix(mvrnorm(n=p, mu=c(0,0), Sigma=cov), p, 2)
      beta = matrix(beta_w[,1],p,1)
      w = matrix(beta_w[,2],p,1)
      delta = rnorm(p, mean=0, sd = sd_delta)
      w = w + delta
      
      y = X%*%beta + pracma::randn(n,1)
      
      ##compute ridge
      eta = lambda * rho * sqrt(a2_t/a2_s)
      beta_t = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% t(X)%*%y
      beta_1D = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% (t(X)%*%y + n*lambda*matrix(w, p, 1))
      beta_2D = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% (t(X)%*%y + n*eta*matrix(w, p, 1))
      
      #inner loop, generate random test data: x_test, y_test
      x_test = pracma::randn(1000,p)%*%S
      y_test = x_test%*%beta + pracma::randn(1000,1)
      
      y_hat_t = x_test%*%beta_t
      y_hat_1D = x_test%*%beta_1D
      y_hat_2D = x_test%*%beta_2D
      mean((y_test - y_hat_t)^2)
      mean((y_test - y_hat_1D)^2)
      mean((y_test - y_hat_2D)^2)
      
      pred_err_t[i,1] = mean((y_test - y_hat_t)^2)
      pred_err_1D[i,1] = mean((y_test - y_hat_1D)^2)
      pred_err_2D[i,1] = mean((y_test - y_hat_2D)^2)
    }
    pred_risk_t[k] = mean(pred_err_t)
    pred_risk_1D[k] = mean(pred_err_1D)
    pred_risk_2D[k] = mean(pred_err_2D)
  }
  return(list(pred_risk_t=pred_risk_t, pred_risk_1D=pred_risk_1D, pred_risk_2D=pred_risk_2D))
}



compute_ST = function(w, t, gamma, grid_size=1e5){
  #w,t - input spectral distribution is a mixture of point masses H = sum delta_{t_i} * w_i
  
  #the v-grid
  v = pracma::linspace(1/grid_size, 1e3, grid_size)
  z = zeros(grid_size,1)
  for(i in 1:grid_size){
    z[i] = -1/v[i] + gamma * sum(w*t/(1 + v[i]*t))
  }
  
  #find the region where z<0, this corresponds to lambda>0
  v = v[z<0]
  z = z[z<0]
  lambda = -z
  
  ind = which((lambda<10) & (lambda>1e-2))
  lambda = lambda[ind]
  v = v[ind]
  z = z[ind]
  m = v/gamma + (1/gamma-1)/z
  
  #compute m',v'
  L = length(lambda)
  v_prime = zeros(L,1)
  for(i in 1:L){
    v_prime[i] = 1 / (1/v[i]^2 - gamma*sum(w*t^2/(1 + t*v[i])^2))
  }
  
  m_prime = v_prime/gamma - (1/gamma-1)/z^2
  
  return(list(lambda=lambda, 
              m=m, 
              v=v, 
              m_prime=m_prime, 
              v_prime=v_prime))
}


#compute risk
#estim_risk - estimation error of ridge on the lambda grid, MSE(beta, beta_est)
#pred_risk - prediction error of ridge on the lambda grid, MSE(y, y_est)
n=10
p=15
rho = 0.6
rate = 1
n_lambda = 50
set.seed(123); sd_delta = sqrt(abs(pracma::randn(p,1)))/2
sigma2 = 1
Vt = 1
Vs = 1
a2_t = p*Vt
a2_s = p*Vs
g = pracma::linspace(1/(2*p), 1-1/(2*p), p)
g = 1/rate*log(1/g)
Sigma = diag(g)
lambda = linspace(0.01, 2.5, n_lambda)^2
t = eig(Sigma)
w = ones(p,1)/p
gamma = p/n

#empirical prediction error
empirical_ridge = simulate_ridge_risk_cov(Sigma, n, p, rho, lambda, sd_delta, Vt=1, Vs=1, sigma2=1)
pred_risk_t = empirical_ridge$pred_risk_t
pred_risk_1D = empirical_ridge$pred_risk_1D
pred_risk_2D = empirical_ridge$pred_risk_2D

plot(sqrt(lambda), pred_risk_t, type = "S", 
     #xlim = c(min(sqrt(lambda)), max(sqrt(lambda))),
     #ylim = c(min(pred_risk_2D, pred_risk_t, na.rm=T), max(pred_risk_2D, pred_risk_t, na.rm=T)),
     xlim = c(0, 2.5), ylim = c(0, 12),
     ylab = 'Prediction error',
     main = 'Correlation=0.6')
lines(sqrt(lambda), pred_risk_1D, type = "S", col='blue')
lines(sqrt(lambda), pred_risk_2D, type = "S", col='red')
legend('bottomright', c('target-only empirical','1D empirical','2D empirical'), col=c('black','blue','red'), lty=c('solid','solid','solid'))

#theoretical prediction error
ST = compute_ST(w, t, gamma)
lambda_th = ST$lambda
m = ST$m
v = ST$v
m_prime = ST$m_prime
v_prime = ST$v_prime
eta = rho * lambda_th * sqrt(a2_t/a2_s)

pred_risk_th_t = (1 + (lambda_th*a2_t/gamma-1)*(1-lambda_th*v_prime/v))/(lambda_th*v)
pred_risk_th_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s - 2*lambda_th*lambda_th*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
pred_risk_th_2D = sigma2 + (lambda_th^2*a2_t + eta^2*a2_s - 2*lambda_th*eta*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)

C_L = min(sd_delta)*p
C_U = max(sd_delta)*p

L_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s + lambda_th^2*C_L - 2*lambda_th*lambda_th*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
U_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s + lambda_th^2*C_U - 2*lambda_th*lambda_th*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)

L_2D = sigma2 + (lambda_th^2*a2_t + eta^2*a2_s + eta^2*C_L - 2*lambda_th*eta*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
U_2D = sigma2 + (lambda_th^2*a2_t + eta^2*a2_s + eta^2*C_U - 2*lambda_th*eta*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)

##match empirical & theoretical
maxi = min(max(lambda), max(lambda_th))
mini = max(min(lambda), min(lambda_th))

ind1 = which( (lambda <=maxi) & (lambda >=mini))
lambda = lambda[ind1]
pred_risk_t = pred_risk_t[ind1]
pred_risk_1D = pred_risk_1D[ind1]
pred_risk_2D = pred_risk_2D[ind1]

ind2 = which((lambda_th <=maxi) & (lambda_th >=mini))
lambda_th = lambda_th[ind2]
L_1D = L_1D[ind2]
L_2D = L_2D[ind2]
U_1D = U_1D[ind2]
U_2D = U_2D[ind2]
pred_risk_th_t = pred_risk_th_t[ind2]
pred_risk_th_1D = pred_risk_th_1D[ind2]
pred_risk_th_2D = pred_risk_th_2D[ind2]


plot(sqrt(lambda), pred_risk_t, type = "S", 
     xlim = c(0, 2.5), ylim = c(5, 20),
     ylab = 'Prediction error', xlab = expression(sqrt(lambda)),
     main = 'Correlation=0.6')
lines(sqrt(lambda), pred_risk_1D, type = "S", col='blue')
#lines(sqrt(lambda_th), L_1D, lty = "dashed", col='blue')
#lines(sqrt(lambda_th), U_1D, lty = "dashed", col='blue')
lines(sqrt(lambda), pred_risk_2D, type = "S", col='red')
#lines(sqrt(lambda_th), L_2D, lty = "dashed", col='red')
#lines(sqrt(lambda_th), U_2D, lty = "dashed", col='red')
polygon(c(sqrt(lambda_th), rev(sqrt(lambda_th))), c(L_2D, rev(U_2D)),
        col=adjustcolor( "red", alpha.f = 0.1), border=NA)
polygon(c(sqrt(lambda_th), rev(sqrt(lambda_th))), c(L_1D, rev(U_1D)),
        col = adjustcolor( "blue", alpha.f = 0.1), border=NA)
legend('bottomright', c('target-only empirical','1D empirical','2D empirical'), col=c('black','blue','red'), lty=c('solid','solid','solid'))






rate = 1
plot_list = list()
order=1
sigma2 = 1
ylabs = c("Prediction Error","Prediction Error", "Prediction Error","Prediction Error", "","", "","", "","", "","")
for(r2 in c(0.3, 0.6, 0.9)){
  for(gamma in c(0.5,2)){
    for(a in c(1,2)){
      if(a==1){
        Vt = 1^2
        Vs = 0.5^2
      }else if(a==2){
        Vt = (1/2)^2
        Vs = (0.9/2)^2
      }
      if(gamma==2){
        if(Vs==(0.9/2)^2){
          plot_list[[order]] = NA
          order = order + 1
          next
        }
        n = 50
        p = 100
        a2_t = p*Vt
        a2_s = p*Vs
        Vratio = Vs/Vt
        yrange = c(5,46)
        xrange = c(0.1,1)
      }else if(gamma==0.5){
        n = 50
        p = 25
        a2_t = p*Vt
        a2_s = p*Vs
        Vratio = Vs/Vt
        if(Vs==(0.5)^2){
          yrange = c(1.5,6)
          xrange = c(0.1,.6)
        }else{
          yrange = c(1,3.1)
        }
        xrange = c(0.1,1)
      }
      #load(paste0("/Users/tiangu/OneDrive - Harvard University/Tian Gu-Shared/ridge/error_band_ration_2_r2_",r2,"_Vratio_", Vratio,".Rdata"))
      load(paste0(import_path, "error_band_gamma_",gamma,"_r2_",r2,"_Vt_", Vt,"_Vs_",Vs,".Rdata"))
      pred_risk_t = tosave$pred_risk_t
      pred_risk_1D = tosave$pred_risk_1D
      pred_risk_2D = tosave$pred_risk_2D
      pred_risk_th_t = tosave$pred_risk_th_t
      pred_risk_th_1D = tosave$pred_risk_th_1D
      pred_risk_th_2D = tosave$pred_risk_th_2D
      lambda = tosave$lambda
      lambda_th = tosave$lambda_th
      U_1D = tosave$U_1D
      U_2D = tosave$U_2D
      L_1D = tosave$L_1D
      L_2D = tosave$L_2D
      
      dat = data.frame(lambda=rep(lambda,3), error=c(pred_risk_t, pred_risk_1D, pred_risk_2D), Method=c(rep('Target only empirical',length(pred_risk_t)),rep('distTL empirical',length(pred_risk_t)),rep('angleTL empirical',length(pred_risk_t))))
      dat$Method = factor(dat$Method, levels = c('Target only empirical','distTL empirical','angleTL empirical'))
      #dat_band = data.frame(lambda=rep(lambda_th,3), lower=c(pred_risk_th_t, L_1D, L_2D), upper=c(pred_risk_th_t, U_1D, U_2D),Color=c(rep('Target only empirical',length(pred_risk_th_t)),rep('distTL empirical',length(pred_risk_th_t)),rep('angleTL empirical',length(pred_risk_th_t))))
      dat_band = data.frame(lambda=rep(lambda_th,3), lower=c(pred_risk_th_t, rep(0,length(lambda_th)), L_2D), upper=c(pred_risk_th_t, rep(0,length(lambda_th)), U_2D),Color=c(rep('Target only empirical',length(pred_risk_th_t)),rep('distTL empirical',length(pred_risk_th_t)),rep('angleTL empirical',length(pred_risk_th_t))))
      dat_band$Color = factor(dat_band$Color, levels = c('Target only empirical','distTL empirical','angleTL empirical'))
      test = merge(x=dat, y=dat_band, by="lambda",all=TRUE)
      
      p = ggplot(test, aes(x=sqrt(lambda),y=error, fill=Method, linetype=Method, color=Method)) + 
        geom_line() + 
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) + 
        xlim(xrange) + theme_bw() + scale_y_continuous(limits = yrange) +
        theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.position = 'none') +
        xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
        ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(sqrt(Vt/Vs))~','~gamma~'='~.(gamma)))+
        #geom_vline(xintercept = sqrt(opt_lam)) +
        scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"))+
        scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'), guide = "legend")+
        scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'))
      
      if(r2==0.9){
        #xrange = c(0.45,1.1)
        p = ggplot(test, aes(x=sqrt(lambda),y=error, fill=Method, linetype=Method, color=Method)) +
          geom_line() +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) +
          xlim(xrange) + theme_bw() + scale_y_continuous(limits = yrange) +
          theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14)) +
          xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
          ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(sqrt(Vt/Vs))~','~gamma~'='~.(gamma)))+
          #geom_vline(xintercept = sqrt(opt_lam)) +
          scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"))+
          scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'))+
          scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'), guide = "legend")
      }
      if(Vs==(0.9/2)^2 & r2!=0.9){
        #xrange = c(0.45,1.1)
        p = ggplot(test, aes(x=sqrt(lambda),y=error, fill=Method, linetype=Method, color=Method)) +
          geom_line() +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) +
          xlim(xrange) + theme_bw() + scale_y_continuous(limits = yrange) +
          theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.position = 'none') +
          xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
          ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=10/9,'~gamma~'='~.(gamma)))+
          scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"))+
          scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'))+
          scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'))
      }else if(Vs==(0.9/2)^2 & r2==0.9){
        p = ggplot(test, aes(x=sqrt(lambda),y=error, fill=Method, linetype=Method, color=Method)) +
          geom_line() +
          #geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) +
          xlim(xrange) + theme_bw() + scale_y_continuous(limits = yrange) +
          theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14)) +
          xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
          ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=10/9,'~gamma~'='~.(gamma)))+
          scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"), guide = "legend")+
          #scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'), na.translate=F)+
          scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'), guide = "legend", na.translate=F)+
          scale_size_manual(values=c(25,25,25))
      }
      plot_list[[order]] = p
      order = order + 1
    }
  }
}
(plot_list[[3]] | plot_list[[7]] | plot_list[[11]])/
  (plot_list[[1]] | plot_list[[5]] | plot_list[[9]])/
  (plot_list[[2]] | plot_list[[6]] | plot_list[[10]])


(plot_list[[1]] | plot_list[[5]] | plot_list[[9]])/
  (plot_list[[2]] | plot_list[[6]] | plot_list[[10]])/
  (plot_list[[4]] | plot_list[[8]] | plot_list[[12]])











# library(ggplot2)
# dat = data.frame(lambda=rep(sqrt(lambda),3), 
#                  error=c(pred_risk_t, pred_risk_1D, pred_risk_2D), 
#                  Method=c(rep('Target only empirical',length(pred_risk_t)),rep('distTL empirical',length(pred_risk_t)),rep('angleTL empirical',length(pred_risk_t))))
# dat$Method = factor(dat$Method, levels = c('Target only empirical','distTL empirical','angleTL empirical'))
# 
# dat_band = data.frame(lambda=rep(sqrt(lambda_th),3), 
#                       lower=c(pred_risk_th_t, L_1D, L_2D), 
#                       upper=c(pred_risk_th_t, U_1D, U_2D),
#                       Color=c(rep('Target only empirical',length(pred_risk_th_t)),rep('distTL empirical',length(pred_risk_th_t)),rep('angleTL empirical',length(pred_risk_th_t))))
# dat_band$Color = factor(dat_band$Color, levels = c('Target only empirical','distTL empirical','angleTL empirical'))
# 
# test = merge(x=dat, y=dat_band, by="lambda",all=TRUE)
# #test$Method = ifelse(test$Color=='Target only empirical', 'Target only theory', test$Method)
# #test$Method = factor(test$Method, levels = c('Target only empirical','Target only theory','distTL empirical','angleTL empirical'))
# 
# ggplot(test, aes(x=lambda,y=error, fill=Method, linetype=Method, color=Method)) + 
#   geom_line() + 
#   geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) + 
#   ylim(c(5, 30)) + xlim(c(0, 2.5)) + theme_bw() +
#   theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14)) +
#   xlab(expression(sqrt(lambda))) + ylab('Prediction Error') +ggtitle('Correlation=0.6')+
#   scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"))+
#   scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'), guide = "legend")+
#   scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'))



