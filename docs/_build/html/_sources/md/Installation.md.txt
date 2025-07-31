# Installation

scATS has been developed with `R 4.0.0` and the following packages are needed to be installed in `R` (R scripts in `scATS` automatically check to see if an R library is already installed and then install those that are needed. So no need for manual preinstallation):

|        module        |    version   |
| -------------------- | ------------ |
| data.table           | >= 1.15.4    |
| fitdistrplus         | >= 1.2-1     |
| GenomicAlignments    | >= 1.22.1    |
| GenomicFeatures      | >= 1.38.2    |
| GenomicRanges        | >= 1.38.0    |
| ggbio                | >= 1.42.0    |
| mclust               | >= 6.1.1     |
| mixtools             | >= 2.0.0     |
| parallel             | >= 3.6.0     |
| rtracklayer          | >= 1.46.0    |
| Seurat               | >= 4.4.0     |
| SummarizedExperiment | >= 1.16.1    |


To install `scATS`, you have two options: either install directly from GitHub or use the compressed source file:
```r
# Install from GitHub if remotes package is not installed
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("LuChenLab/r-scATS/scATS/")
```

Alternatively, you can install `scATS` using the source file downloaded from the [repository](https://github.com/LuChenLab/r-scATS) :
```r
# Install scATS from a downloaded source file
R CMD INSTALL scATS_0.5.4.tar.gz
```
