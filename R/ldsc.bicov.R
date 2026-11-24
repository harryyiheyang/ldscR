#' Bivariate LD Score Regression (LDSC) for Cross-Trait Genetic Covariance
#'
#' `ldsc.bicov` estimates the genetic covariance between two traits by weighted
#' least squares of \eqn{z_{1j} z_{2j}} on \eqn{(1, \ell_j \sqrt{N_{1j} N_{2j}} / M)}.
#'
#' @param gwas1,gwas2 `data.frame`s with columns `SNP`, `Zscore`, `N`, aligned
#'   to the same alleles.
#' @param h21,h22 Heritability estimates of the two traits (used in the weights
#'   and to report `rg`).
#' @param LDSC LD-score panel as in [ldsc.univ()] (columns `SNP`, `LDSC`,
#'   optional `wLD`); all panel SNPs.
#' @param sampling.time,nblock,block.counts Block bootstrap, as in [ldsc.univ()].
#' @param M Number of panel SNPs (`NULL`: rows of `LDSC`).
#' @param two.step Chi-square cutoff: the intercept is estimated from SNPs with
#'   both \eqn{\chi^2_1, \chi^2_2 <} `two.step`, the slope from all SNPs.
#'   `NULL` = one step.
#' @param z.max Truncation of both Z-scores; `NULL` = \eqn{\sqrt{\max(80, 0.001 N_{max})}}.
#'
#' @details For bivariate normal \eqn{(z_1, z_2)} with variances
#' \eqn{v_1, v_2} and covariance \eqn{\mu},
#' \eqn{\mathrm{Var}(z_1 z_2) = v_1 v_2 + \mu^2}. The weights are therefore
#' \deqn{w_j^{-1} = \max(1, \ell^{w}_j)\,[(1 + N_{1j} h_{21} \ell_j / M)(1 + N_{2j} h_{22} \ell_j / M) + \mu_j^2],}
#' with \eqn{\mu_j = \mathrm{ecov} + \mathrm{gcov}\, \ell_j \sqrt{N_{1j}N_{2j}}/M}
#' updated twice; the univariate intercepts in \eqn{v_1, v_2} are set to 1.
#'
#' @return A `data.frame` with `ecov`, `ecov.se`, `gcov`, `gcov.se`, `rg`
#'   (`gcov / sqrt(h21 h22)`, no SE: use [ldsc.rg()]), `M`, `M.reg`.
#' @export
ldsc.bicov = function(gwas1, gwas2, h21, h22, LDSC, sampling.time = 0, nblock = 500, M = NULL, two.step = 30, z.max = NULL, block.counts = NULL) {
  gwas1 <- as.data.frame(gwas1); gwas2 <- as.data.frame(gwas2); LDSC <- as.data.frame(LDSC)
  if (is.null(M)) M <- length(unique(LDSC$SNP))
  if (!"wLD" %in% colnames(LDSC)) LDSC$wLD <- LDSC$LDSC
  common <- LDSC$SNP[LDSC$SNP %in% gwas1$SNP & LDSC$SNP %in% gwas2$SNP]
  LDSC <- LDSC[match(common, LDSC$SNP), ]
  gwas1 <- gwas1[match(common, gwas1$SNP), ]
  gwas2 <- gwas2[match(common, gwas2$SNP), ]
  M.reg <- length(common)
  if (is.null(z.max)) z.max <- .zmax_default(c(gwas1$N, gwas2$N))

  z1 <- .truncate_z(gwas1$Zscore, z.max); z2 <- .truncate_z(gwas2$Zscore, z.max)
  y <- z1 * z2
  l <- LDSC$LDSC * sqrt(gwas1$N) * sqrt(gwas2$N) / M
  X <- cbind(1, l)
  lw <- 1 / pmax(1, LDSC$wLD)
  base_var <- (1 + LDSC$LDSC * gwas1$N / M * max(h21, 0)) * (1 + LDSC$LDSC * gwas2$N / M * max(h22, 0))
  use2 <- !is.null(two.step)
  S1 <- if (use2) z1^2 < two.step & z2^2 < two.step else rep(TRUE, M.reg)
  all <- rep(TRUE, M.reg)
  .xp <- function(w, rows) list(XtX = crossprod(X[rows, , drop = FALSE], X[rows, , drop = FALSE] * w[rows]),
                                Xty = as.vector(crossprod(X[rows, , drop = FALSE], y[rows] * w[rows])))

  w <- lw / base_var
  b <- .fit_ab(.xp(w, S1), .xp(w, all), use2)
  for (i in 1:2) {
    mu <- b[1] + b[2] * l
    w <- lw / (base_var + mu^2)
    b <- .fit_ab(.xp(w, S1), .xp(w, all), use2)
  }
  ecov <- b[1]; gcov <- b[2]

  boot <- sampling.time > 0 || !is.null(block.counts)
  ecov.se <- gcov.se <- NA_real_
  if (boot) {
    IndexList <- .block_index(M.reg, nblock)
    b1 <- .block_xprod(X, y, w * S1, IndexList)
    b2 <- .block_xprod(X, y, w, IndexList)
    cnt <- .boot_counts(length(IndexList), sampling.time, block.counts)
    est <- t(apply(cnt, 1, function(cc) .fit_ab(.sum_xprod(b1, cc), .sum_xprod(b2, cc), use2)))
    ecov.se <- sd(est[, 1]); gcov.se <- sd(est[, 2])
  }
  rg <- if (h21 > 0 && h22 > 0) gcov / sqrt(h21 * h22) else NA_real_
  result <- data.frame(ecov = ecov, ecov.se = ecov.se, gcov = gcov, gcov.se = gcov.se, rg = rg, M = M, M.reg = M.reg)
  if (boot) result <- .new_ldsc_bootstrap_result(result, list(ecov = est[, 1], gcov = est[, 2], block.counts = cnt))
  result
}

#' Genetic correlation by LDSC with a joint block bootstrap
#'
#' Fits [ldsc.univ()] for both traits and [ldsc.bicov()] on the SNPs common to
#' both GWAS and the panel, reusing the same bootstrap draws in all three fits,
#' so that \eqn{r_g = \mathrm{gcov}/\sqrt{h^2_1 h^2_2}} gets a bootstrap SE.
#'
#' @inheritParams ldsc.bicov
#' @return A `data.frame` with `h21`, `h21.se`, `h22`, `h22.se`, `gcov`,
#'   `gcov.se`, `ecov`, `ecov.se`, `rg`, `rg.se`, `M`, `M.reg`.
#' @export
ldsc.rg = function(gwas1, gwas2, LDSC, sampling.time = 500, nblock = 500, M = NULL, two.step = 30, z.max = NULL) {
  gwas1 <- as.data.frame(gwas1); gwas2 <- as.data.frame(gwas2); LDSC <- as.data.frame(LDSC)
  if (is.null(M)) M <- length(unique(LDSC$SNP))
  common <- LDSC$SNP[LDSC$SNP %in% gwas1$SNP & LDSC$SNP %in% gwas2$SNP]
  gwas1 <- gwas1[match(common, gwas1$SNP), ]
  gwas2 <- gwas2[match(common, gwas2$SNP), ]
  if (is.null(z.max)) z.max <- .zmax_default(c(gwas1$N, gwas2$N))
  B <- length(.block_index(length(common), nblock))
  cnt <- .boot_counts(B, sampling.time)
  u1 <- ldsc.univ(gwas1, LDSC, nblock = nblock, M = M, two.step = two.step, z.max = z.max, block.counts = cnt)
  u2 <- ldsc.univ(gwas2, LDSC, nblock = nblock, M = M, two.step = two.step, z.max = z.max, block.counts = cnt)
  bc <- ldsc.bicov(gwas1, gwas2, u1$h2, u2$h2, LDSC, nblock = nblock, M = M, two.step = two.step, z.max = z.max, block.counts = cnt)
  rg.b <- bc$bootstrap$gcov / sqrt(pmax(u1$bootstrap$h2, 0) * pmax(u2$bootstrap$h2, 0))
  rg.b[!is.finite(rg.b)] <- NA
  data.frame(h21 = u1$h2, h21.se = u1$h2.se, h22 = u2$h2, h22.se = u2$h2.se,
             gcov = bc$gcov, gcov.se = bc$gcov.se, ecov = bc$ecov, ecov.se = bc$ecov.se,
             rg = bc$rg, rg.se = sd(rg.b, na.rm = TRUE), M = M, M.reg = length(common))
}
