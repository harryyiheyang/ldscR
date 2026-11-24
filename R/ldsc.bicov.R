#' Bivariate LD Score Regression (LDSC): genetic covariance and correlation
#'
#' `ldsc.bicov` regresses \eqn{z_{1j} z_{2j}} on \eqn{(1, \ell_j \sqrt{N_{1j} N_{2j}} / M)}
#' by weighted least squares. It always fits the two univariate regressions
#' ([ldsc.univ()]) on the same SNPs inside the same bootstrap: every replicate
#' re-estimates \eqn{h^2_1}, \eqn{h^2_2} and the covariance, so
#' \eqn{r_g = \mathrm{gcov}/\sqrt{h^2_1 h^2_2}} and all SEs come from one set of draws. The computation is in C++ (`src/ldsc_univ_bicov.cpp`).
#'
#' @param gwas1,gwas2 `data.frame`s with columns `SNP`, `Zscore`, `N`, aligned
#'   to the same alleles.
#' @param LDSC LD-score panel as in [ldsc.univ()] (columns `SNP`, `LDSC`,
#'   optional `wLD`); all panel SNPs.
#' @param sampling.time,nblock Block bootstrap, as in [ldsc.univ()].
#' @param M Number of panel SNPs (`NULL`: rows of `LDSC`).
#' @param two.step Chi-square cutoff: the intercepts are estimated from SNPs
#'   with \eqn{\chi^2 <} `two.step` (both traits for the cross-trait fit), the
#'   slopes from all SNPs. `NULL` = one step.
#' @param z.max Truncation of both Z-scores; `NULL` = \eqn{\sqrt{\max(80, 0.001 N_{max})}}.
#'
#' @details For bivariate normal \eqn{(z_1, z_2)} with variances
#' \eqn{v_1, v_2} and covariance \eqn{\mu},
#' \eqn{\mathrm{Var}(z_1 z_2) = v_1 v_2 + \mu^2}. The weights are therefore
#' \deqn{w_j^{-1} = \max(1, \ell^{w}_j)\,[(1 + N_{1j} h_{21} \ell_j / M)(1 + N_{2j} h_{22} \ell_j / M) + \mu_j^2],}
#' with \eqn{\mu_j = \mathrm{ecov} + \mathrm{gcov}\, \ell_j \sqrt{N_{1j}N_{2j}}/M}
#' updated twice; the univariate intercepts in \eqn{v_1, v_2} are set to 1.
#'
#' @return A `data.frame` with `ecov`, `ecov.se`, `gcov`, `gcov.se`, `h21`,
#'   `h21.se`, `h22`, `h22.se`, `rg`, `rg.se`, `M`. `rg` is set to 0 when
#'   \eqn{h^2_1 \le 0} or \eqn{h^2_2 \le 0} (no genetic variance), in the full fit and in
#'   every bootstrap replicate. With bootstrap,
#'   `$bootstrap` holds the replicate vectors and `block.counts`.
#' @export
ldsc.bicov = function(gwas1, gwas2, LDSC, sampling.time = 0, nblock = 500, M = NULL, two.step = 30, z.max = NULL) {
  gwas1 <- as.data.frame(gwas1); gwas2 <- as.data.frame(gwas2); LDSC <- as.data.frame(LDSC)
  if (is.null(M)) M <- length(unique(LDSC$SNP))
  if (!"wLD" %in% colnames(LDSC)) LDSC$wLD <- LDSC$LDSC
  common <- LDSC$SNP[LDSC$SNP %in% gwas1$SNP & LDSC$SNP %in% gwas2$SNP]
  LDSC <- LDSC[match(common, LDSC$SNP), ]
  gwas1 <- gwas1[match(common, gwas1$SNP), ]
  gwas2 <- gwas2[match(common, gwas2$SNP), ]
  if (is.null(z.max)) z.max <- .zmax_default(c(gwas1$N, gwas2$N))
  cnt <- .boot_counts(.n_groups(length(common), nblock), sampling.time)
  r <- ldsc_bicov_cpp(as.numeric(gwas1$Zscore), as.numeric(gwas2$Zscore), as.numeric(gwas1$N), as.numeric(gwas2$N),
                      as.numeric(LDSC$LDSC), as.numeric(LDSC$wLD), as.numeric(M),
                      z.max, .two_step_arg(two.step), as.integer(nblock), cnt, 2L, 1e-12)
  .rg <- function(g, a, b) ifelse(a > 0 & b > 0, g / sqrt(pmax(a, 0) * pmax(b, 0)), 0)
  rg <- .rg(r$b[2], r$h21, r$h22)
  se <- .ldsc_bootstrap_sd(r$boot)
  h.se <- c(NA_real_, NA_real_); rg.se <- NA_real_; rg.b <- NULL
  if (sampling.time > 0) {
    h.se <- c(stats::sd(r$boot1[, 2]), stats::sd(r$boot2[, 2]))
    rg.b <- .rg(r$boot[, 2], r$boot1[, 2], r$boot2[, 2])
    rg.se <- stats::sd(rg.b)
  }
  result <- data.frame(ecov = r$b[1], ecov.se = se[1], gcov = r$b[2], gcov.se = se[2],
                       h21 = r$h21, h21.se = h.se[1], h22 = r$h22, h22.se = h.se[2],
                       rg = rg, rg.se = rg.se, M = M)
  if (sampling.time > 0) {
    boot <- list(ecov = r$boot[, 1], gcov = r$boot[, 2], h21 = r$boot1[, 2], h22 = r$boot2[, 2], rg = rg.b, block.counts = cnt)
    result <- .new_ldsc_bootstrap_result(result, boot)
  }
  result
}
