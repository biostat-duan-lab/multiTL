## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true"),
  knitr::opts_chunk$set(purl = NOT_CRAN)
)

## ----setup--------------------------------------------------------------------
library(multiTL)
library(glmnet)
library(stats)
library(viRandomForests)
library(randomForest)

## -----------------------------------------------------------------------------
### Prepare data
load(paste0("../data/data_commute.rda"))
X = data_commute$X
y = data_commute$y
X.src = data_commute$X.src
y.src = data_commute$y.src
n.src = data_commute$n.src
w.src = data_commute$w.src
r = data_commute$r
S = data_commute$S

### Run COMMUTE
res = commute (X, y, X.src, y.src, n.src, w.src, r, S)

## -----------------------------------------------------------------------------
### Prepare data
load(paste0("../data/data_transRF.rda"))
X = data_transRF$X
y = data_transRF$y
rf.src = data_transRF$rf.src
S = data_transRF$S

### Run transRF
res = transRF (X, y,X.test=NULL, rf.src, S)

## -----------------------------------------------------------------------------
### Prepare data
load(paste0("../data/data_angleTL.rda"))
X = data_angleTL$X
y = data_angleTL$y
w.src = data_angleTL$w.src

### Run angleTL
res_angle = angleTL (X, y, w.src)

### Run distTL
res_dist = distTL (X, y, w.src)

