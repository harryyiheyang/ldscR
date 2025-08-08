#' Bivariate LD Score Regression (LDSC) for Cross‑Trait Genetic Covariance
#'
#' `ldsc.bicov` estimates the genetic covariance between two traits using
#' LD Score Regression (LDSC). The function harmonizes two GWAS
#' summary-statistic tables with an LD‑score table, fits a weighted least
#' squares (WLS) regression of \eqn{z_{1j} z_{2j}} on the LD regressor, and
#' performs one (optional) reweighting step for efficiency. Block bootstrap
#' can be used to obtain standard errors.
#'
#' @param gwas1 A `data.frame`/`data.table` of GWAS summary statistics for trait 1.
#'   Required columns:
#'   \itemize{
#'     \item `SNP` — SNP identifier
#'     \item `Zscore` — Z‑score of the SNP–trait association
#'     \item `N` — per‑SNP sample size
#'   }
#' @param gwas2 A `data.frame`/`data.table` of GWAS summary statistics for trait 2,
#'   with the same required columns as `gwas1`.
#' @param LDSC A `data.frame`/`data.table` containing LD scores. Required columns:
#'   \itemize{
#'     \item `SNP` — SNP identifier
#'     \item `LDSC` — LD score
#'   }
#' @param h21 Numeric. Heritability estimate of trait 1 (used in weights).
#' @param h22 Numeric. Heritability estimate of trait 2 (used in weights).
#' @param nblock Integer. Number of blocks for block bootstrap.
#' @param sampling.time Integer. Number of bootstrap replicates. If `0`, only
#'   point estimates are returned.
#'
#' @details
#' Let \eqn{M} be the number of SNPs after harmonization. The regressor is
#' \eqn{\ell_j = \mathrm{LDSC}_j \sqrt{N_{1j}} \sqrt{N_{2j}} / M}. The initial
#' WLS uses weights
#' \deqn{w_j^{-1} = \bigl(1 + \mathrm{LDSC}_j N_{1j} h_{21} / M \bigr)
#'                     \bigl(1 + \mathrm{LDSC}_j N_{2j} h_{22} / M \bigr),}
#' followed by an efficiency reweighting with the mean term
#' \eqn{\mu_j = \mathrm{ecov} + \mathrm{gcov} \cdot \ell_j} via
#' \eqn{w_j^{-1} \leftarrow w_j^{-1} + 2 \mu_j^2}.
#'
#' @return A `data.frame` with:
#' \itemize{
#'   \item `ecov` — estimated intercept (environmental covariance)
#'   \item `ecov.se` — standard error of `ecov` (bootstrap; `NA` if no bootstrap)
#'   \item `gcov` — estimated genetic covariance (slope)
#'   \item `gcov.se` — standard error of `gcov` (bootstrap; `NA` if no bootstrap)
#'   \item `M` — number of SNPs analyzed
#' }
#'
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom data.table setDT setkey
#' @export
ldsc.bicov = function(gwas1, gwas2, h21, h22, LDSC, nblock = 500, sampling.time = 0) {
  data.table::setDT(LDSC); data.table::setkey(LDSC, SNP)

  # --------------------------- Harmonize inputs ---------------------------
  t0 <- Sys.time()
  SNPintersect <- intersect(gwas1$SNP, gwas2$SNP)
  gwas1 <- gwas1[gwas1$SNP %in% SNPintersect, ]
  gwas2 <- gwas2[gwas2$SNP %in% SNPintersect, ]
  gwas1 <- LDSC[gwas1, nomatch = 0]
  gwas2 <- LDSC[gwas2, nomatch = 0]
  M <- nrow(gwas1)
  t0 <- Sys.time() - t0
  print("Processing data"); print(t0)

  # --------------------------- Helper: one WLS fit ------------------------
  .wls_fit <- function(l, z, w) {
    X <- cbind(1, l)
    XtX <- CppMatrix::matrixMultiply(t(X), X * w)
    Xty <- CppMatrix::matrixVectorMultiply(t(X), z * w)
    as.vector(solve(XtX) %*% Xty)
  }

  # --------------------------- Initial estimator --------------------------
  t1 <- Sys.time()
  z <- gwas1$Zscore * gwas2$Zscore
  l <- gwas1$LDSC * sqrt(gwas1$N) * sqrt(gwas2$N) / M
  w <- 1 / (1 + gwas1$LDSC * gwas1$N / M * h21) /
    (1 + gwas2$LDSC * gwas2$N / M * h22)
  beta <- .wls_fit(l, z, w)
  gcov.ini <- beta[2]; ecov.ini <- beta[1]
  t1 <- Sys.time() - t1
  print("Initial Genetic Covariance Estimate"); print(t1)

  # --------------------------- Reweighting step ---------------------------
  t2 <- Sys.time()
  mu <- l * gcov.ini + ecov.ini
  w  <- 1 / ((1 + gwas1$LDSC * gwas1$N / M * h21) /
               (1 + gwas2$LDSC * gwas2$N / M * h22) + 2 * mu^2)
  beta <- .wls_fit(l, z, w)
  gcov <- beta[2]; ecov <- beta[1]

  mu <- l * gcov + ecov
  w  <- 1 / ((1 + gwas1$LDSC * gwas1$N / M * h21) /
               (1 + gwas2$LDSC * gwas2$N / M * h22) + 2 * mu^2)
  beta <- .wls_fit(l, z, w)
  gcov <- beta[2]; ecov <- beta[1]
  t2 <- Sys.time() - t2
  print("Genetic Covariance Reweighting"); print(t2)

  # --------------------------- Block bootstrap ---------------------------
  if (sampling.time > 0) {
    t3 <- Sys.time()
    blocksize <- floor(M / nblock)
    IndexMatrix <- matrix(seq_len(nblock * blocksize), nrow = blocksize, ncol = nblock)
    gap <- M - (nblock * blocksize)
    h2.vec <- intercept.vec <- numeric(sampling.time)

    for (i in seq_len(sampling.time)) {
      indcol <- sample.int(ncol(IndexMatrix), ncol(IndexMatrix), replace = TRUE)
      ind <- as.vector(IndexMatrix[, indcol]) + ifelse(gap > 1, floor(runif(1, 0, gap)), 0)

      z1 <- z[ind]; l1 <- l[ind]; w1 <- w[ind]
      b  <- .wls_fit(l1, z1, w1)
      h2.vec[i] <- b[2]
      intercept.vec[i] <- b[1]
    }
    t3 <- Sys.time() - t3
    print("Standard Error Estimation"); print(t3)

    ecov.se <- sqrt(mean((intercept.vec - ecov)^2))
    gcov.se <- sqrt(mean((h2.vec - gcov)^2))
  } else {
    gcov.se <- ecov.se <- NA_real_
  }

  data.frame(ecov = ecov, ecov.se = ecov.se,
             gcov = gcov, gcov.se = gcov.se,
             M = M)
}
