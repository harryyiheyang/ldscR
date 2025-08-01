% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldsc.bicov.R
\name{ldsc.bicov}
\alias{ldsc.bicov}
\title{Single-Variate Linkage Disequilibrium Score Regression (LDSC)}
\usage{
ldsc.bicov(gwas1, gwas2, h21, h22, LDSC, nblock = 500, sampling.time = 0)
}
\arguments{
\item{gwas1}{A data.frame containing gwas summary statistics for the first trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N).}

\item{gwas2}{A data.frame containing gwas summary statistics for the second trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N).}

\item{LDSC}{A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.}

\item{nblock}{The number of blocks for bootstrap-based standard error estimation.}

\item{sampling.time}{The number of block bootstrap. If sampling.time=0, only genetic covariance estimate will be returned.}
}
\value{
A data.frame containing genetic covariance estimates and their standard errors, along with the estimation error covariance estimate and its standard error.
}
\description{
The `ldsc.bicov` function performs bi-variate genetic covariance using Linkage Disequilibrium Score Regression (LDSC) analysis. It is designed to estimate heritability from gwas1 summary statistics, accounting for linkage disequilibrium (LD) between SNPs. The function harmonizes gwas1 data with LD scores and applies non-linear optimization to estimate heritability.
}
\details{
The `ldsc.univ` function is particularly useful for single-trait heritability analysis using LDSC. It requires well-formatted gwas1 summary statistics and an understanding of LDSC methodology. The function includes an initial estimation phase followed by a reweighting step for improved efficiency and accuracy.
}
