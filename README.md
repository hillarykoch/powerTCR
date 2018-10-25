# powerTCR: modeling the clone size distribution of the TCR repertoire

This is an R package for fitting the discrete gamma-GPD spliced threshold model to a distribution of clone sizes. The package contains tools needed to perform all of the analyses found in __our paper (pending)__. 

## Installation

Install and load this package from BioConductor by typing in R:

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("powerTCR")
```

or take it directly from from GitHub with:

```{r}
library(devtools)
install_github("hillarykoch/powerTCR")
library(powerTCR)
```

## Getting going

See the [package vignette](/vignettes/powerTCR.Rmd) for a detailed walkthrough of package features.

## Citation

If you found the powerTCR package useful, please cite our paper:

Koch, H., Starenki, D., Cooper, S.J., Myers, R.M., Li, Q. (2018) powerTCR: a model-based approach to comparative analysis of the clone size distribution of the T cell receptor repertoire, PLoS Computational Biology.

While our paper is yet to appear in print, you may see our [preprint on bioRvix:](https://www.biorxiv.org/content/early/2018/04/07/297119)

Citation for the powerTCR package is:

Koch H (2018). powerTCR: Model-Based Comparative Analysis of the TCR Repertoire. R package version 1.1.4.


