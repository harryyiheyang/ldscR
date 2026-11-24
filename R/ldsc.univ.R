#' Univariate LD Score Regression (LDSC) for SNP-Heritability
#'
#' `ldsc.univ` estimates SNP-heritability for a single trait by weighted least
#' squares of \eqn{Z_j^2} on \eqn{(1, \ell_j N_j / M)}, with an optional block
#' bootstrap for standard errors.
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
#' @param block.counts Optional bootstrap count matrix (replicates x blocks) to
#'   reuse the draws of another fit (see [ldsc.rg()]).
#'
#' @details Weights: \eqn{w_j = 1/\max(1, \ell^{w}_j) / (1 + \hat h^2 N_j \ell_j / M)^2},
#' started from \eqn{\hat h^2 = 0.1} and updated three times; the intercept in
#' the weights is set to 1.
#'
#' @return A `data.frame` with `intercept`, `intercept.se`, `h2`, `h2.se`, `M`
#'   (panel SNPs) and `M.reg` (regression SNPs). With bootstrap, `$bootstrap`
#'   holds the replicate vectors and `block.counts`.
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom data.table setDT setkey
#' @export
ldsc.univ = function(gwas, LDSC, sampling.time = 0, nblock = 500, M = NULL, two.step = 30, z.max = NULL, block.counts = NULL) {
  gwas <- as.data.frame(gwas); LDSC <- as.data.frame(LDSC)
  if (is.null(M)) M <- length(unique(LDSC$SNP))
  if (!"wLD" %in% colnames(LDSC)) LDSC$wLD <- LDSC$LDSC
  common <- gwas$SNP[gwas$SNP %in% LDSC$SNP]
  gwas <- gwas[match(common, gwas$SNP), ]
  LDSC <- LDSC[match(common, LDSC$SNP), ]
  M.reg <- nrow(gwas)
  if (is.null(z.max)) z.max <- .zmax_default(gwas$N)

  y <- .truncate_z(gwas$Zscore, z.max)^2
  l <- LDSC$LDSC * gwas$N / M
  X <- cbind(1, l)
  lw <- 1 / pmax(1, LDSC$wLD)
  use2 <- !is.null(two.step)
  S1 <- if (use2) y < two.step else rep(TRUE, M.reg)
  all <- rep(TRUE, M.reg)
  .xp <- function(w, rows) list(XtX = crossprod(X[rows, , drop = FALSE], X[rows, , drop = FALSE] * w[rows]),
                                Xty = as.vector(crossprod(X[rows, , drop = FALSE], y[rows] * w[rows])))

  w <- lw / (1 + 0.1 * l)^2
  beta <- .fit_ab(.xp(w, S1), .xp(w, all), use2)
  for (i in 1:3) {
    w <- lw / (1 + l * max(beta[2], 0))^2
    beta <- .fit_ab(.xp(w, S1), .xp(w, all), use2)
  }
  intercept <- beta[1]; h2 <- beta[2]

  boot <- sampling.time > 0 || !is.null(block.counts)
  h2.se <- intercept.se <- NA_real_
  if (boot) {
    IndexList <- .block_index(M.reg, nblock)
    b1 <- .block_xprod(X, y, w * S1, IndexList)
    b2 <- .block_xprod(X, y, w, IndexList)
    cnt <- .boot_counts(length(IndexList), sampling.time, block.counts)
    est <- t(apply(cnt, 1, function(cc) .fit_ab(.sum_xprod(b1, cc), .sum_xprod(b2, cc), use2)))
    intercept.se <- sd(est[, 1]); h2.se <- sd(est[, 2])
  }
  result <- data.frame(intercept = intercept, intercept.se = intercept.se, h2 = h2, h2.se = h2.se, M = M, M.reg = M.reg)
  if (boot) result <- .new_ldsc_bootstrap_result(result, list(intercept = est[, 1], h2 = est[, 2], block.counts = cnt))
  result
}
