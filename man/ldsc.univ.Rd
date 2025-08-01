% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldsc.univ.R
\name{ldsc.univ}
\alias{ldsc.univ}
\title{Single-Variate Linkage Disequilibrium Score Regression (LDSC)}
\usage{
ldsc.univ(gwas, LDSC, nblock = 500, sampling.time = 500)
}
\arguments{
\item{gwas}{A data.frame containing GWAS summary statistics for a single trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).}

\item{LDSC}{A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.}

\item{nblock}{The number of blocks for bootstrap-based standard error estimation.}

\item{sampling.time}{The number of block bootstrap. If sampling.time=0, only genetic covariance estimate will be returned.}
}
\value{
A data.frame containing heritability estimates and their standard errors, along with the intercept and its standard error.
}
\description{
The `ldsc.uni` function performs single-variate Linkage Disequilibrium Score Regression (LDSC) analysis. It is designed to estimate heritability from GWAS summary statistics, accounting for linkage disequilibrium (LD) between SNPs. The function harmonizes GWAS data with LD scores and applies non-linear optimization to estimate heritability.
}
\details{
The `ldsc.univ` function is particularly useful for single-trait heritability analysis using LDSC. It requires well-formatted GWAS summary statistics and an understanding of LDSC methodology. The function includes an initial estimation phase followed by a reweighting step for improved efficiency and accuracy.
}
