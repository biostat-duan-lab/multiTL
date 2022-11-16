#' @title area_under_curve
#' @description Calculating AUC
#' @param y,y.prob estimated y and actual y
#' @return AUC
#' @importFrom pROC roc

area_under_curve <- function(y, y.prob){
  pROC::auc(pROC::roc(y, as.numeric(y.prob)))
}

