
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Phase2.2CATE

<!-- badges: start -->
<!-- badges: end -->

The goal of Phase2.2CATE is to perform Federated Adaptive Causal
Estimation of Conditional Vaccine Effect.

## Installation

This package should be run locally (i.e. NOT on docker). You can install
the development version of Phase2.2CATE like so:

``` r
remotes::install_github('clara-lea/Phase2.2CATE')
```

## Running the analysis

These are the lines of code to run CATE

``` r
library(Phase2.2CATE)
input.path = "" #specify the path to the folder containing 2.2 data
output.path = "" #specify where you would like the results to be saved
siteid = "" #specify siteid

run_analysis(input.path, output.path, siteid)
```

After running the analysis, please send the result file directly to
Clara-Lea.
