#' Simulated data for COMMUTE
#' @name data_commute
#' @docType data
NULL

# library(pROC)
#
# #simulation SNP data K population
# Coef.gen<- function(s, h, K, sig.delta, p, exact=TRUE, intercept){
#
#   beta <- c(seq(0.4, 0.1, length.out=(s/2)), -seq(0.4, 0.1, length.out=(s/2)), rep(0, p - s))
#   w <- matrix(rep(beta, K), p, K)   ###h^*=0
#   if(exact==TRUE){
#     samp <- sample(1:p, h, replace=F)
#     sign <- sample(c(-1,1), h, replace=T)
#     for(k in 1:K){
#       w[samp, k] <- w[samp, k] + sign*rep(as.numeric(sig.delta[k]), h)
#     }
#   }else{
#     total_p = 1:p
#     samp0 <- sample(total_p, min(h), replace=F)
#     sign <- sample(c(-1,1), h, replace=T)
#     samp <- samp0
#     w[samp, 1] <- w[samp, 1] + sign*rep(as.numeric(sig.delta), h[1])
#     for(k in 2:K){
#       if(h[k] > h[k-1]){
#         h_diff = h[k] - h[k-1]
#         samp_more <- sample(total_p[-samp], h_diff, replace=F)
#         sign_more <- sample(c(-1,1), h_diff, replace=T)
#         w[samp, k] <- w[samp, k-1]
#         w[samp_more, k] <- w[samp_more, k-1] + sign_more*rep(as.numeric(sig.delta), h_diff)
#       }else{
#         w[, k] <- w[, k-1]
#       }
#     }
#   }
#
#   beta = c(intercept[1], beta)
#   w = rbind(intercept[-1], w)
#
#   return(list(beta=beta, w=w))
# }
#
#
# library(sim1000G)
#
# ##### generate X pool
# vcf_file = file.path(system.file("examples", package = "sim1000G"),"region.vcf.gz")
# vcf = readVCF(vcf_file, maxNumberOfVariants = 500, min_maf = 0.1, max_maf = 0.9)
# downloadGeneticMap( 4 )
# readGeneticMap( chromosome = 4 )
# startSimulation(vcf)
# Sig1 = SIM$haplodata$cor
# X.pool = NULL
# for(row_rep in 1:10){
#   X.pool.tmp = NULL
#   for(col_rep in 1:4){
#     print(col_rep)
#     SIM$reset(); id = c()
#     for(i in 1:2000) id[i] = SIM$addUnrelatedIndividual()
#     genotypes = SIM$gt1[id,] + SIM$gt2[id,]
#     X.pool.tmp = cbind(X.pool.tmp, genotypes)
#     col_rep = col_rep + 1
#   }
#   X.pool = rbind(X.pool, X.pool.tmp)
#   row_rep = row_rep + 1
# }
#
# exact = FALSE             #TRUE=S1; FALSE=S2
# K = 3                     #total number of source population
# p = 2000
# s = 100
# n0 = 100                  #target population size
# nt = c(1000, 1500, 2000)    #source population size
# # nt = c(2000, 3000, 4000)    #source population size
# n.test = 1000
# intercept = rep(1, K+1)
#
# n.tar <- n0
# n.src <- nt
# sig.delta=0.8
# h = c(30,30,30)
#
# coef.all <- Coef.gen(s=s, h=h, K=K, sig.delta, p=p, exact, intercept)
# beta.true <- as.numeric(coef.all$beta)
# B=cbind(beta.true, coef.all$w)
#
# sample.tar = sample(1:nrow(X.pool), n.tar, replace = FALSE)
# X.tar = X.pool[sample.tar,]
# y.tar<- rbinom(n.tar, size=1, prob=logistic(B[1,1]+X.tar%*%B[-1, 1]))
# X.pool = X.pool[-sample.tar,]
# X.src <- y.src <- list()
# for(k in 1:K){
#   sample.k = sample(1:nrow(X.pool), n.src[k], replace = FALSE)
#   X.src[[k]] = X.pool[sample(1:nrow(X.pool), n.src[k], replace = FALSE),]
#   y.src[[k]]<- rbinom(n.src[k], size=1, prob=logistic(B[1, k+1]+X.src[[k]]%*%B[-1, k+1]))
#   X.pool = X.pool[-sample.k,]
# }
#
# data_commute = list(X=data_commute$X,y=data_commute$y,X.src=data_commute$X.src,y.src=data_commute$y.src,n.src=data_commute$n.src,w.src=data_commute$w.src,r=25,S=2)
# X.src = X.src[1:2]
# y.src = y.src[1:2]
# w.src = list()
# w.src[[1]] = coef.all$w[,3]
# n.src=nt
#
# # commute(X.tar, y.tar, X.src[1:2], y.src[1:2], n.src=nt, w.src=coef.all$beta[3], r=1, S=2)
