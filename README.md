# powerTCR: modeling the clone size distribution of the TCR repertoire

This is an R package for fitting the discrete gamma-GPD spliced threshold model to a distribution of clone sizes. The package contains tools needed to perform all of the analyses found in __our paper (pending)__. 

## Installation

Install and load this package from BioConductor by typing in R:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("powerTCR")
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

Paper yet to appear.

See our [preprint on bioRvix:](https://www.biorxiv.org/content/early/2018/04/07/297119)
Hillary Koch, Dmytro Starenki, Sara J Cooper, Richard M Myers, Qunhua Li
bioRxiv 297119; doi: https://doi.org/10.1101/297119


