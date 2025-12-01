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

## Step 1: Load LD Score References

This package supports two types of LD Score reference panels:

1.  HapMap3-based references (~1.1-1.7M SNPs)

- Traditional sparse reference panels
- Built-in data: Hapmap3_EURLDSC, Hapmap3_EASLDSC, etc.
- Suitable for standard heritability analyses
- Faster computation due to fewer SNPs

2.  High-density 7M references (~7M SNPs)

- Based on UK Biobank imputed data via SBayesRC
- Provides enhanced accuracy for complex trait analyses
- Available for single ancestries (EUR, EAS, AFR) and cross-population
  combinations
- Requires download due to large file sizes

### Option A: Use Built-in HapMap3 References

``` r
library(ldscR)
# Load built-in European HapMap3 LD scores
data("Hapmap3")
data("Hapmap3_EURLDSC")
head(Hapmap3)
```

    ##           SNP     A1     A2
    ##        <char> <char> <char>
    ## 1:  rs2185539      T      C
    ## 2:  rs3131972      G      A
    ## 3: rs12184325      T      C
    ## 4:  rs3131969      G      A
    ## 5:  rs3131967      C      T
    ## 6:  rs3131962      G      A

``` r
head(Hapmap3_EURLDSC)
```

    ##                   SNP     LDSC
    ## rs3131972   rs3131972 2.520038
    ## rs3131969   rs3131969 2.523700
    ## rs1048488   rs1048488 2.522668
    ## rs12562034 rs12562034 1.028920
    ## rs4040617   rs4040617 2.495103
    ## rs4970383   rs4970383 2.020633

### Option B: Download High-Density 7M References

``` r
library(ldscR)
# Download 7M European LD scores (recommended for enhanced accuracy)
download_7MLDSC("EUR")
```

    ## Downloading 7M LD Score reference for EUR ancestry...
    ##   Downloading 7m_eurldsc.parquet ...
    ##   Downloading 7m_snp.parquet ...
    ##   Loading SNP information...

    ## Warning: Potentially unsafe or invalid elements have been discarded from R metadata.
    ## ℹ Type: "externalptr"
    ## → If you trust the source, you can set `options(arrow.unsafe_metadata = TRUE)` to preserve them.

    ##   Loading LD scores...
    ##   Loaded 7356518 SNPs with LD scores as SevenMillon_EURLDSC 
    ##   SNP information available as SevenMillon_SNP with 7356518 SNPs
    ## Download complete! Data loaded to workspace.

``` r
# This loads SevenMillon_EURLDSC and SevenMillon_SNP into your workspace

head(SevenMillon_EURLDSC)
```

    ##           SNP     LDSC
    ## 1  rs12132974 8.193143
    ## 2  rs12134490 8.193950
    ## 3  rs17276806 8.195587
    ## 4 rs139867617 8.191437
    ## 5   rs7526310 9.609900
    ## 6  rs72631880 8.605243

``` r
head(SevenMillon_SNP)
```

    ##            SNP   CHR     BP     A1     A2     Freq Block MarkerName
    ##         <char> <int>  <int> <char> <char>    <num> <int>     <char>
    ## 1:  rs12132974     1 801661      T      C 0.079275     1   1:801661
    ## 2:  rs12134490     1 801680      C      A 0.079225     1   1:801680
    ## 3:  rs17276806     1 801858      T      C 0.079300     1   1:801858
    ## 4: rs139867617     1 802856      T      C 0.079275     1   1:802856
    ## 5:   rs7526310     1 804759      T      C 0.123250     1   1:804759
    ## 6:  rs72631880     1 805556      A      T 0.085425     1   1:805556

``` r
# Alternative: Download cross-population references
# download_7MLDSC("EURxEAS")  # European x East Asian
# download_7MLDSC("EURxAFR")  # European x African
```

## Step 2: Download Example Data

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
```

## Step 3: Align SNPs

``` r
gwaslist <- filter_align(
  gwas_data_list = gwaslist,
  ref_panel = Hapmap3[, c("SNP", "A1", "A2")]
)
```

    ## Adjusting effect allele according to reference panel...
    ## Finding common SNPs...
    ## Aligning data to common SNPs and ordering...
    ## Filtering complete.

## Step 4: Estimate SNP Heritability (Univariate LDSC)

``` r
fit1 <- ldsc.univ(
  gwas = gwaslist$driving,
  LDSC = Hapmap3_EURLDSC,
  nblock = 500,
  sampling.time = 200
)
```

    ## [1] "Processing data"
    ## Time difference of 0.4858479 secs
    ## [1] "Initial Heritability Estimate"
    ## Time difference of 0.03685904 secs
    ## [1] "Heritability Estimation (Reweighted)"
    ## Time difference of 0.04275489 secs
    ## [1] "Standard Error Estimation"
    ## Time difference of 13.38615 secs

``` r
fit2 <- ldsc.univ(
  gwas = gwaslist$computer,
  LDSC = Hapmap3_EURLDSC,
  nblock = 500,
  sampling.time = 200
)
```

    ## [1] "Processing data"
    ## Time difference of 0.6447089 secs
    ## [1] "Initial Heritability Estimate"
    ## Time difference of 0.026932 secs
    ## [1] "Heritability Estimation (Reweighted)"
    ## Time difference of 0.411679 secs
    ## [1] "Standard Error Estimation"
    ## Time difference of 13.14892 secs

``` r
fit1  # Heritability and intercept for driving time
```

    ##   intercept intercept.se         h2        h2.se       M
    ## 1  1.024993  0.004432506 0.03044821 0.0006913167 1038356

``` r
fit2  # Heritability and intercept for computer use time
```

    ##   intercept intercept.se         h2        h2.se       M
    ## 1   1.04081  0.005311228 0.07046025 0.0007884596 1038356

## Step 5: Estimate Genetic Covariance (Bivariate LDSC)

``` r
fit3 <- ldsc.bicov(
  gwas1 = gwaslist$driving,
  gwas2 = gwaslist$computer,
  LDSC = Hapmap3_EURLDSC,
  nblock = 500,
  sampling.time = 200,
  h21 = fit1$h2,
  h22 = fit2$h2
)
```

    ## [1] "Processing data"
    ## Time difference of 1.244697 secs
    ## [1] "Initial Genetic Covariance Estimate"
    ## Time difference of 0.02895212 secs
    ## [1] "Genetic Covariance Reweighting"
    ## Time difference of 0.04821205 secs
    ## [1] "Standard Error Estimation"
    ## Time difference of 14.46051 secs

``` r
fit3  # Genetic and sample-overlap covariance between traits
```

    ##         ecov     ecov.se        gcov      gcov.se       M
    ## 1 0.05856721 0.005963673 0.000789279 0.0008040988 1038356

## Step 6: Construct Covariance Matrices

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
