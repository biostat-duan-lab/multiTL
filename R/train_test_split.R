#' @title train_test_split
#' @description Spliting the training dataset and the testing dataset
#' @param x x
#' @param y y
#' @param size The ratio to split for training data
#' @return A list of x_train, y_train, x_test, y_test

train_test_split <- function(x, y, size){
  n = length(y)
  train = sample(1:n, size*n)
  test = 1:n
  test = test[!test %in% train]
  x_train = x[train,]
  y_train = y[train]
  x_test = x[test,]
  y_test = y[test]
  return(list(x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test))
}
