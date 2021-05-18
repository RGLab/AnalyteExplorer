# AnalyteExplorer

<!-- badges: start -->
[![R-CMD-check](https://github.com/RGLab/AnalyteExplorer/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/AnalyteExplorer/actions)
<!-- badges: end -->

The goal of AnalyteExplorer is to pre-process data for the `AnalyteExplorer` module in [ImmuneSpace](https://immunespace.org/)

## Installation

You can install the development version of AnalyteExplorer from [GitHub](https://CRAN.R-project.org) with:

``` r
# install.packages("remotes")
remotes::install_github("RGLab/AnalyteExplorer")
```

## Example

``` r
library(AnalyteExplorer)
labkey.url.base <- "https://www.immunespace.org"
labkey.url.path <- "/AnalyteExplorer"
data <- process_data("gene_expression")
check_data(data)
res <- update_table("gene_expression", data)
```
