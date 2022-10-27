#' @title transRF
#' @description a transfer learning random forest framework
#' @param X Variables from the target. The variables need to be completely the same set and in the same order as variables used in the source model.
#' @param y Response from the target.
#' @param X.test If no testing data is availble, automatically split 20\% target data for testing.
#' @param rf.src Trained random forest model from the source.
#' @param S Feature importance score from the source model.
#' @param p.val The percent of spliting training data for validation. The default value is 10\%.
#' @return A list of predicted y and y.test if X.test is not provided.
#' @importFrom viRandomForests viRandomForests
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export


transRF = function(X, y, X.test=NULL, rf.src, S, p.val=0.1){

  ####### If no testing data is availble, split 20% target data for testing
  if(is.null(X.test)){
    sample0 <- sample(c(TRUE, FALSE), nrow(X), replace=TRUE, prob=c(0.8, 0.2))
    X.tar0 <- X[sample0, ]
    X.test <- X[!sample0, ]
    y.tar0 <- y[sample0]
    y.test <- y[!sample0]
  }else{
    X.tar0 = X
    y.tar0 = y
    y.test <- rep(NA, dim(X.test)[1])
  }
  dat.test <- data.frame(y.test=y.test, X.test)

  var_list = attr(rf.src$terms,"term.labels")
  num = length(var_list)

  ###### Spliting p.val% (default 10%) training data fot validation
  sample1 <- sample(c(TRUE, FALSE), nrow(X.tar0), replace=TRUE, prob=c(1-p.val, p.val))
  X.tar <- X.tar0[sample1, ]
  X.val <- X.tar0[!sample1, ]
  y.tar <- y.tar0[sample1]
  y.val <- y.tar0[!sample1]

  dat.tar = data.frame(y.tar=y.tar, X.tar)
  dat.val = data.frame(y.val=y.val, X.val)


  ###### Target-only: fit RandomForest
  colnames(dat.tar)[1]='y.tar'
  rf.m0 <- randomForest(y.tar ~ ., data=dat.tar, ntree=500)
  colnames(dat.test)=colnames(dat.tar)
  y.m0.pred = predict(rf.m0, dat.test)

  ###### Source model prediction
  colnames(dat.test)[1:num+1]=var_list
  y.src.pred = predict(rf.src , dat.test)

  ###### Model 1: fit viRandomForest y.tar ~ X.tar with source.score
  rf.m1 <- viRandomForests(y.tar ~ ., data=dat.tar, ntree=500, fprob=S, keep.forest=TRUE, importance=TRUE)
  colnames(dat.test)=colnames(dat.tar)
  y.m1.pred = predict(rf.m1, dat.test)

  ###### Model 2: fit viRandomForest y.delta ~ X.tar
  colnames(dat.tar)[1:num+1]=var_list
  y.src.hat = predict(rf.src , dat.tar)
  y.delta = y.tar - y.src.hat

  dat.train2 = cbind(y.delta, dat.tar[,-1])
  colnames(dat.train2)[1]='y.delta'
  rf.delta <- viRandomForests(y.delta ~ ., data=dat.train2, ntree=500, fprob=NULL, keep.forest=TRUE, importance=TRUE)
  colnames(dat.train2)[1:num+1]=var_list
  colnames(dat.test)[1:num+1]=var_list
  y.delta.pred = predict(rf.delta, dat.test)
  y.m2.pred = y.src.pred + y.delta.pred

  ###### Model 3: fit y.tar ~ X.tar + source.pred
  dat.y.src.train = cbind(dat.tar, y.src.hat)
  dat.y.src.test = data.frame(dat.test, y.src.pred)
  colnames(dat.y.src.test) = colnames(dat.y.src.train)
  rf.m3 <- viRandomForests(y.tar ~ ., data=dat.y.src.train, ntree=500,  fprob=c(rep(1,length(S)),2), keep.forest=TRUE, importance=TRUE)
  y.m3.pred = predict(rf.m3, dat.y.src.test)

  ###### TransRF: emsemble of target-only and Models 1-3
  colnames(dat.val)=colnames(dat.tar)
  y0 = predict(rf.m0, dat.val)
  y1 = predict(rf.m1, dat.val)

  colnames(dat.val)[1:num+1]=var_list
  y2 = predict(rf.delta, dat.val) + predict(rf.src, dat.val)
  y3 = predict(rf.m3, data.frame(dat.val, y.src.hat=predict(rf.src , dat.val)))
  y4 = predict(rf.src , dat.val)

  weight.ensemble = coef(lm(y.val~y0+y1+y2+y3))
  y.transrf = cbind(1,y.m0.pred,y.m1.pred, y.m2.pred, y.m3.pred) %*% weight.ensemble

  output = list(y.transrf = y.transrf, y.test=y.test) #if user provides X.test, then y.test=NULL
  return(output)
}




