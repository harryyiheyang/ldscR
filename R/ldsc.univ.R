#' Univariate LD Score Regression (LDSC) for SNP-Heritability
#'
#' `ldsc.univ` estimates SNP-heritability for a single trait by weighted least
#' squares of \eqn{Z_j^2} on \eqn{(1, \ell_j N_j / M)}, with an optional block
#' bootstrap for standard errors. The computation is in C++ (`src/ldsc_univ_bicov.cpp`).
#'
#' @param gwas A `data.frame` with columns `SNP`, `Zscore`, `N` (alleles must
#'   already be aligned to the LD panel, e.g. with [filter_align()]).
#' @param LDSC A `data.frame` with columns `SNP` and `LDSC` (LD score). An
#'   optional column `wLD` gives the LD score used in the over-counting weight
#'   \eqn{1/\max(1, \ell^{w}_j)} (default: `LDSC`). This must be the whole
#'   LD panel: SNPs missing from the GWAS stay in the table.
#' @param sampling.time Number of bootstrap replicates (`0` = point estimates only).
#' @param nblock Number of contiguous SNP blocks for the bootstrap.
#' @param M Number of SNPs over which the LD scores were summed. `NULL` uses
#'   the number of SNPs in `LDSC` (the panel), not the number of regression
#'   SNPs: \eqn{E[\chi^2_j] = a + N h^2 \ell_j / M_{panel}} holds whichever
#'   SNPs the GWAS reports.
#' @param two.step Chi-square cutoff for the two-step estimator: the intercept
#'   is estimated from SNPs with \eqn{\chi^2 <} `two.step`, then the slope from
#'   all SNPs with the intercept fixed. `NULL` fits both in one step.
#' @param z.max Z-scores are truncated (not removed) to `[-z.max, z.max]`.
#'   `NULL` uses \eqn{\sqrt{\max(80, 0.001 N_{max})}}, the ldsc chi-square cap.
#'
#' @details Weights: \eqn{w_j = 1/\max(1, \ell^{w}_j) / (1 + \hat h^2 N_j \ell_j / M)^2},
#' started from \eqn{\hat h^2 = 0.1} and updated three times; the intercept in
#' the weights is set to 1.
#'
#' @return A `data.frame` with `intercept`, `intercept.se`, `h2`, `h2.se`, `M`
#'   (panel SNPs) and `M.reg` (regression SNPs). With bootstrap, `$bootstrap`
#'   holds the replicate vectors and `block.counts`.
#' @useDynLib ldscR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
ldsc.univ = function(gwas, LDSC, sampling.time = 0, nblock = 500, M = NULL, two.step = 30, z.max = NULL) {
  gwas <- as.data.frame(gwas); LDSC <- as.data.frame(LDSC)
  if (is.null(M)) M <- length(unique(LDSC$SNP))
  if (!"wLD" %in% colnames(LDSC)) LDSC$wLD <- LDSC$LDSC
  common <- gwas$SNP[gwas$SNP %in% LDSC$SNP]
  gwas <- gwas[match(common, gwas$SNP), ]
  LDSC <- LDSC[match(common, LDSC$SNP), ]
  if (is.null(z.max)) z.max <- .zmax_default(gwas$N)
  cnt <- .boot_counts(.n_groups(nrow(gwas), nblock), sampling.time)
  r <- ldsc_univ_cpp(as.numeric(gwas$Zscore), as.numeric(gwas$N), as.numeric(LDSC$LDSC), as.numeric(LDSC$wLD),
                     as.numeric(M), z.max, .two_step_arg(two.step), as.integer(nblock), cnt, 3L, 1e-12)
  se <- .ldsc_bootstrap_sd(r$boot)
  result <- data.frame(intercept = r$b[1], intercept.se = se[1], h2 = r$b[2], h2.se = se[2], M = M, M.reg = nrow(gwas))
  if (sampling.time > 0) result <- .new_ldsc_bootstrap_result(result, list(intercept = r$boot[, 1], h2 = r$boot[, 2], block.counts = cnt))
  result
}
