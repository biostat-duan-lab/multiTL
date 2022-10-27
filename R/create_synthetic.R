#' @title create_synthetic
#' @description Create synthetic data under federated settings
#' @param K Number of source sites that individual data cannot be shared
#' @param X.tar X from target data
#' @param n.src Source sample size
#' @param r A positive integer
#' @param B A matrix of beta and w, beta is in the 1st column
#' @importFrom stats rbinom
#' @return A list of synthetic X and Y
#' @export
#'
create_synthetic <- function(K, X.tar, n.src, r, B){
  n.syn = round(n.src*r)     #total number of synthetic data for each source: r*source sample

  X.syn <- y.syn <- list()
  for(k in 1:K){
    print(k)

    X.syn.k = X.tar[sample(1:nrow(X.tar), size = n.syn[k], replace=TRUE),]
    y.syn.k = rbinom(n.syn[k], 1, prob = logistic(c(cbind(1,X.syn.k) %*% B[[k]])))

    X.syn[[k]] = X.syn.k
    y.syn[[k]] = y.syn.k
  }
  return(list(X.syn=X.syn, y.syn=y.syn))
}
