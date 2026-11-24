#' Trans-Ancestry Bivariate Linkage Disequilibrium Score Regression (LDSC)
#'
#' The `ldsc.trans` function performs bivariate genetic covariance estimation
#' using Linkage Disequilibrium Score Regression (LDSC) in a trans-ancestry setting.
#' It estimates the genetic covariance between two traits from different ancestries,
#' using a cross-ancestry LD score for the regression regressor and ancestry-specific
#' LD scores for the variance weights. The function also supports block bootstrap
#' for standard error estimation.
#'
#' @param gwas1 A data.frame or data.table containing GWAS summary statistics
#'   for the first trait (ancestry 1). Must include:
#'   \itemize{
#'     \item \code{SNP}: SNP identifier (string)
#'     \item \code{Zscore}: Z-score for the SNP-trait association
#'     \item \code{N}: Sample size for the SNP
#'   }
#' @param gwas2 A data.frame or data.table containing GWAS summary statistics
#'   for the second trait (ancestry 2), with the same column requirements as \code{gwas1}.
#' @param h21 Numeric. Heritability estimate of the first trait.
#' @param h22 Numeric. Heritability estimate of the second trait.
#' @param LDSC1 A data.frame or data.table containing ancestry-1 LD scores for
#'   variance weighting. Must include:
#'   \itemize{
#'     \item \code{SNP}: SNP identifier
#'     \item \code{LDSC}: LD score value for ancestry 1
#'   }
#' @param LDSC2 A data.frame or data.table containing ancestry-2 LD scores for
#'   variance weighting. Must include:
#'   \itemize{
#'     \item \code{SNP}: SNP identifier
#'     \item \code{LDSC}: LD score value for ancestry 2
#'   }
#' @param LDSC_Tran A data.frame or data.table containing cross-ancestry LD scores
#'   for the regression regressor. Must include:
#'   \itemize{
#'     \item \code{SNP}: SNP identifier
#'     \item \code{LDSC}: Cross-ancestry LD score
#'   }
#' @param nblock Integer. Number of blocks used for block bootstrap standard error
#'   estimation.
#' @param sampling.time Integer. Number of bootstrap replicates. If set to 0,
#'   only point estimates will be returned.
#'
#' @param M Number of panel SNPs of the cross-ancestry LD score (`NULL`: rows of `LDSC_Tran`).
#'
#' @details Weights \eqn{w_j^{-1} = (1 + \ell_{1j} N_{1j} h_{21}/M)(1 + \ell_{2j} N_{2j} h_{22}/M) + \mu_j^2},
#'   \eqn{\mu_j} the fitted mean, updated twice. Computed in C++ (`src/ldsc_univ_bicov.cpp`).
#' @return A data.frame with `ecov`, `ecov.se`, `gcov`, `gcov.se`, `M` (panel SNPs).
#' @export
ldsc.trans = function(gwas1, gwas2, h21, h22, LDSC1, LDSC2, LDSC_Tran, sampling.time = 0, nblock = 500, M = NULL) {
  gwas1 <- as.data.frame(gwas1); gwas2 <- as.data.frame(gwas2)
  LDSC1 <- as.data.frame(LDSC1); LDSC2 <- as.data.frame(LDSC2); LDSC_Tran <- as.data.frame(LDSC_Tran)
  if (is.null(M)) M <- length(unique(LDSC_Tran$SNP))
  common <- LDSC_Tran$SNP[LDSC_Tran$SNP %in% gwas1$SNP & LDSC_Tran$SNP %in% gwas2$SNP &
                          LDSC_Tran$SNP %in% LDSC1$SNP & LDSC_Tran$SNP %in% LDSC2$SNP]
  g1 <- gwas1[match(common, gwas1$SNP), ]; g2 <- gwas2[match(common, gwas2$SNP), ]
  cnt <- .boot_counts(.n_groups(length(common), nblock), sampling.time)
  r <- ldsc_trans_cpp(as.numeric(g1$Zscore), as.numeric(g2$Zscore), as.numeric(g1$N), as.numeric(g2$N),
                      as.numeric(LDSC1$LDSC[match(common, LDSC1$SNP)]), as.numeric(LDSC2$LDSC[match(common, LDSC2$SNP)]),
                      as.numeric(LDSC_Tran$LDSC[match(common, LDSC_Tran$SNP)]), as.numeric(M), h21, h22,
                      as.integer(nblock), cnt, 2L, 1e-12)
  se <- .ldsc_bootstrap_sd(r$boot)
  data.frame(ecov = r$b[1], ecov.se = se[1], gcov = r$b[2], gcov.se = se[2], M = M)
}
