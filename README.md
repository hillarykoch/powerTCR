# powerTCR: modeling the clone size distribution of the TCR repertoire

This is an R package for fitting the discrete gamma-GPD spliced threshold model to a distribution of clone sizes. The package contains tools needed to perform all of the analyses found in our paper, *powerTCR: a model-based approach to comparative analysis of the clone size distribution of the T cell receptor repertoire.*

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

You may read our methods paper at [PLoS Computational Biology:](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006571)

If you found the powerTCR package useful, please cite our paper:

Koch H, Starenki D, Cooper SJ, Myers RM, Li Q (2018) powerTCR: A model-based approach to comparative analysis of the clone size distribution of the T cell receptor repertoire. PLOS Computational Biology 14(11): e1006571. https://doi.org/10.1371/journal.pcbi.1006571

Citation for the powerTCR package is:

Koch H (2018). powerTCR: Model-Based Comparative Analysis of the TCR Repertoire. R package version 1.1.4.


