Estimating Heritability and Genetic Covariance via LDSC
================
Yihe Yang (<yxy1234@case.edu>)

## Overview

This tutorial demonstrates how to perform single- and bivariate Linkage
Disequilibrium Score Regression (LDSC) to estimate:

- SNP heritability of two behavioral traits: **driving time** and
  **computer use time**;
- Genetic and environmental covariance between them.

The example uses preprocessed GWAS summary statistics, originally
derived from UK Biobank questionnaires, and LD scores based on European
HapMap3 reference.

------------------------------------------------------------------------

## Installation and Setup

``` r
devtools::install_github("harryyiheyang/ldscR")
```

## Step 1: Load Example Data

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ldscR)
data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode = "wb")
gwaslist <- readRDS(temp_file)
unlink(temp_file)

# Load reference alleles from HapMap3
data("hapmap3")
```

## Step 2: Align SNPs

``` r
gwaslist <- filter_align(
  gwas_data_list = gwaslist,
  ref_panel = hapmap3[, c("SNP", "A1", "A2")]
)
```

    ## Adjusting effect allele according to reference panel...
    ## Finding common SNPs...
    ## Aligning data to common SNPs and ordering...
    ## Filtering complete.

## Step 3: Estimate SNP Heritability (Univariate LDSC)

``` r
fit1 <- ldsc.univ(
  gwas = gwaslist$driving,
  LDSC = EURLDSC,
  nblock = 500,
  sampling.time = 200
)
```

    ## [1] "Processing data"
    ## Time difference of 1.296875 secs
    ## [1] "Initial Genetic Covariance Estimate"
    ## Time difference of 0.05574298 secs
    ## [1] "Heritability Estimation"
    ## Time difference of 0.1449771 secs
    ## [1] "Standard Error Estimation"
    ## Time difference of 8.578486 secs

``` r
fit2 <- ldsc.univ(
  gwas = gwaslist$computer,
  LDSC = EURLDSC,
  nblock = 500,
  sampling.time = 200
)
```

    ## [1] "Processing data"
    ## Time difference of 0.667803 secs
    ## [1] "Initial Genetic Covariance Estimate"
    ## Time difference of 0.02304792 secs
    ## [1] "Heritability Estimation"
    ## Time difference of 0.09671807 secs
    ## [1] "Standard Error Estimation"
    ## Time difference of 8.654569 secs

``` r
fit1  # Heritability and intercept for driving time
```

    ##   intercept intercept.se         h2       h2.se       M
    ## 1  1.024993   0.01537545 0.03044821 0.002714929 1038356

``` r
fit2  # Heritability and intercept for computer use time
```

    ##   intercept intercept.se         h2       h2.se       M
    ## 1   1.04081   0.01617612 0.07046025 0.003542177 1038356

## Step 4: Estimate Genetic Covariance (Bivariate LDSC)

``` r
fit3 <- ldsc.bicov(
  gwas1 = gwaslist$driving,
  gwas2 = gwaslist$computer,
  LDSC = EURLDSC,
  nblock = 500,
  sampling.time = 200,
  h21 = fit1$h2,
  h22 = fit2$h2
)
```

    ## [1] "Processing data"
    ## Time difference of 0.5761371 secs
    ## [1] "Initial Genetic Covariance Estimate"
    ## Time difference of 0.1039009 secs
    ## [1] "Heritability Estimation"
    ## Time difference of 0.09886599 secs
    ## [1] "Standard Error Estimation"
    ## Time difference of 8.219073 secs

``` r
fit3  # Genetic and sample-overlap covariance between traits
```

    ##         ecov     ecov.se        gcov      gcov.se       M
    ## 1 0.05856721 0.005404381 0.000789279 0.0007017296 1038356

## Step 5: Construct Covariance Matrices

``` r
gcov_matrix <- matrix(
  c(fit1$h2, fit3$gcov,
    fit3$gcov, fit2$h2),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Driving", "Computer"), c("Driving", "Computer"))
)

ecov_matrix <- matrix(
  c(fit1$intercept, fit3$ecov,
    fit3$ecov, fit2$intercept),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Driving", "Computer"), c("Driving", "Computer"))
)

gcov_matrix
```

    ##              Driving    Computer
    ## Driving  0.030448210 0.000789279
    ## Computer 0.000789279 0.070460248

``` r
ecov_matrix
```

    ##             Driving   Computer
    ## Driving  1.02499318 0.05856721
    ## Computer 0.05856721 1.04081016
