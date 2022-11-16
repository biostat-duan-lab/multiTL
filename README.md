# multiTL

The goal of multiTL is to incorporate pre-trained models from source into transfer learning framework. We offer three major methods: commute, transRF and angleTL. For details about the methodology, check out our paper: 
- [1] Gu T, Lee PH, Duan R. “COMMUTE: Communication-efficient transfer learning approach for multi-site risk prediction.” (2022) [MedRxiv](https://www.medrxiv.org/content/10.1101/2022.03.23.22272834v1).
- [2] Gu T, Han Y, Duan R. “A transfer learning approach based on random forest with application to breast cancer prediction in underrepresented
    populations.”  (2022) [Proceedings of Pacific Symposium on Biocomputing (PSB) 2023](https://psb.stanford.edu/callfor/papers/psb23_papers_allv2.pdf).
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

If you are on Linux or you would like to compile from source, you can download the source codes [multiTL_0.1.0.tar.gz](https://github.com/biostat-duan-lab/multiTL/releases/multiTL_0.1.0.tar.gz). Mac users should refer to [this page](https://cran.r-project.org/bin/macosx/tools/) for the various dependencies required. Install then via: 
```r
install.packages("/path/to/downloaded_source.tar.gz", repos=NULL, type="source")
```

If you have `devtools`, you can also type: 
```r
install_github("biostat-duan-lab/multiTL")
```
