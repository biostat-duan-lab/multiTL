#' @title commute
#' @description Derive calibration term estimator and an estimator for the target parameter in each source sites
#' @param X Variables from the target. The variables need to be completely the same set and in the same order as variables used in the model parameter estimators.
#' @param y Response from the target.
#' @param X.src A list of shareable sources X
#' @param y.src A list of shareable sources y
#' @param n.src Sample size of K sources
#' @param w.src A list of source estimators from which individual level data cannot be shared
#' @param r A positive integer
#' @param S The first S commute estimators to be incorporated in the aggregation
#' @importFrom glmnet cv.glmnet
#' @return A list of two beta estimator using inverse weight and logistic regression respectively
#' @export

commute = function(X, y, X.src, y.src, n.src, w.src, r, S){

  n.tar = nrow(X)
  p = ncol(X)
  M = length(X.src)
  K = M + length(w.src)

  if(length(n.src)!=K) stop('The length of n.src does not match.')
  # if(p!=p.src) stop('The dimensions of source and target do not match.')

  ### split training and validation set
  xy_split = train_test_split(X, y, size=0.8)
  x_train = xy_split$x_train
  y_train = xy_split$y_train
  x_val = xy_split$x_test
  y_val = xy_split$y_test


  ### directly use glmnet to fit lasso
  beta.tar <- as.numeric(ST_init(x_train, y_train)$beta0)
  w = list()
  TL = list()
  for(k in 1:M){
    w[[k]] <- as.numeric(ST_init(X.src[[k]], y.src[[k]])$beta0)
    TL[[k]] <- TL_init(x_train, y_train, X.src=X.src[[k]], y.src=y.src[[k]], w=w[[k]])
  }

  for(k in 1:(K-M)){
    w[[k+M]] <- w.src[[k]]
    TL[[k+M]] <- TL_init(x_train, y_train, X.src=NULL, y.src=NULL, w=w[[k+M]])
  }


  delta.TL = list() ###calculate delta for each source
  w.thres = list() #source-only LASSO + threshold
  for(k in 1:K){
    delta.TL[[k]] = thres(TL[[k]]$delta0, n.tar, p) ###add threshold
    w.thres[[k]] = thres(w[[k]], sqrt(n.src[k]), p) ###add threshold
  }

  ### SURE Screening
  auc.tar = c()
  for(k in 1:K){
    pred.y = logistic(x_train%*%w.thres[[k]][-1])
    auc.tmp = area_under_curve(y_train, pred.y)
    auc.tar = c(auc.tar, auc.tmp)
  }

  rank.auc.tmp = c()
  rank.auc = c()
  for (i in 1:K) {
    rank.auc.tmp = c(rank.auc.tmp,which(auc.tar==sort(auc.tar, decreasing = T)[i]))
  }
  # rank.auc = rank.auc.tmp[1:S]

  ### singleTL
  beta.TL.single = list()
  for(k in 1:K){
    beta.TL.single[[k]] = TL[[k]]$beta0 ###add threshold
  }

  ### COMMUTE
  data.syn <- create_synthetic(K, x_train, n.src, r, B=beta.TL.single)
  beta.commute = list()

  for (i in 1:S) {
    rank.auc = rank.auc.tmp[1:i]

    if(any(rank.auc<=M)==FALSE){
      X.tar.s <- do.call(rbind, data.syn$X.syn[rank.auc])
      y.tar.s <- do.call(c, data.syn$y.syn[rank.auc])
      X.tar.adj <- rbind(x_train, X.tar.s)
      y.tar.adj <- rbind(y_train, y.tar.s)
      beta.rank = ST_init(X.tar.adj, y.tar.adj)$beta0

    }else if(any(rank.auc>M)==FALSE){
      rank.auc.src = rank.auc[rank.auc < M+1]
      X.src.s <- X.src[rank.auc.src]
      y.src.s <- y.src[rank.auc.src]
      beta.rank = Trans_global(x_train, y_train, X.src=X.src.s, y.src=y.src.s, delta=delta.TL[rank.auc.src])

    }else{
      rank.auc.tar = rank.auc[rank.auc > M]
      rank.auc.src = rank.auc[rank.auc < M+1]

      X.src.s <- X.src[rank.auc.src]
      y.src.s <- y.src[rank.auc.src]
      X.tar.s <- do.call(rbind, data.syn$X.syn[rank.auc.tar])
      y.tar.s <- do.call(c, data.syn$y.syn[rank.auc.tar])
      X.tar.adj <- rbind(x_train, X.tar.s)
      y.tar.adj <- c(y_train, y.tar.s)

      beta.rank = Trans_global(X.tar.adj, y.tar.adj, X.src=X.src.s, y.src=y.src.s, delta=delta.TL[rank.auc.src])
    }
    beta.commute[[i]]=beta.rank
  }

  ### aggregation
  B.singleTL = do.call(cbind, beta.TL.single)
  B.commute = do.call(cbind, beta.commute)
  B.syn = cbind(B.singleTL, beta.tar, B.commute)

  wt.syn.agg1 <- Agg_fun1(B.syn, x_val, y_val, const=1)
  wt.syn.agg2 <- Agg_fun2(B.syn, x_val, y_val, const=1)
  beta.inv <- B.syn%*%wt.syn.agg1
  beta.logit <- B.syn%*%wt.syn.agg2

  return(list(beta.inv=beta.inv, beta.logit=beta.logit))
}




