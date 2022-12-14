# multiTL

In multiTL, we offer algorithms of three methods to achieve the same goal-- improving the risk prediction performance of a target population with limited samples by incorporating pre-trained models from external source data through transfer learning approaches. Details of these three methods--COMMUTE [1], transRF [2] and angleTL [3], can be found in our preprints and paper:

- [1] Gu T, Lee PH, Duan R. “COMMUTE: Communication-efficient transfer learning approach for multi-site risk prediction.” (2022) [Journal of Biomedical Informatics](https://doi.org/10.1016/j.jbi.2022.104243).
- [2] Gu T, Han Y, Duan R. “A transfer learning approach based on random forest with application to breast cancer prediction in underrepresented
    populations.”  (2022) [Proceedings of Pacific Symposium on Biocomputing (PSB) 2023](https://www.worldscientific.com/doi/pdf/10.1142/9789811270611_0018).
- [3] Gu T, Han Y, Duan R. “Robust angle-based transfer learning in high dimensions”  (2022) [ArXiv](http://arxiv.org/abs/2210.12759).

## Installation

`multiTL` requires the following R packages: 'corpcor', 'MASS', 'stats', 'viRandomForests', 'randomForest', 'pROC', 'glmnet'. Install them by: 

```r
install.packages(c('corpcor', 'MASS', 'stats', 'randomForest', 'pROC', 'glmnet'), dependencies=TRUE)
```

For R package 'viRandomForests', you can download the source codes [viRandomForests_1.0.tar.gz](https://github.com/biostat-duan-lab/multiTL/blob/master/viRandomForests_1.0.tar.gz)

```r
install.packages('/path/to/viRandomForests_1.0.tar.gz', type='source', repo=NULL)
```

If you are on Linux or you would like to compile from source, you can download the source codes [multiTL_0.1.0.tar.gz](https://github.com/biostat-duan-lab/multiTL/blob/master/releases/multiTL_0.1.0.tar.gz). Mac users should refer to [this page](https://cran.r-project.org/bin/macosx/tools/) for the various dependencies required. Install then via: 
```r
install.packages("/path/to/downloaded_source.tar.gz", repos=NULL, type="source")
```
