
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
devtools::install_github("claralea/Phase2.2CATE")
```

## Vaccine data

The vaccine data should be saved in the same folder as 2.2 data with the
file name ‘LocalPatientVaccine.csv’. Its columns are the following:

-   siteid

-   cohort

-   patient_num

-   vaccine_type (type of the first vaccine: Moderna, Pfizer, JJ, Other,
    None)

-   vaccine_date (date of the first vaccine: if no vaccine indicate
    1900-01-01)

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
