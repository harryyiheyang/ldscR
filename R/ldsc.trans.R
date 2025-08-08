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
#' @return A data.frame with the following columns:
#'   \itemize{
#'     \item \code{ecov}: Estimated intercept (environmental covariance)
#'     \item \code{ecov.se}: Standard error of the intercept
#'     \item \code{gcov}: Estimated genetic covariance
#'     \item \code{gcov.se}: Standard error of the genetic covariance
#'     \item \code{M}: Number of SNPs used in the analysis
#'   }
#'
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom data.table setDT setkey
#'
#' @export

ldsc.trans = function(gwas1, gwas2, h21, h22, LDSC1, LDSC2, LDSC_Tran,
                      nblock = 500, sampling.time = 0) {
  # Expect columns:
  # gwas1/gwas2: SNP, Zscore, N
  # LDSC_Tran:   SNP, LDSC12              (cross-ancestry LD score for the regressor)
  # LDSC1:       SNP, LDSC1               (ancestry-1 within LD score, for weights)
  # LDSC2:       SNP, LDSC2               (ancestry-2 within LD score, for weights)

  #----------------------- Harmonize once on a single skeleton -----------------------
  data.table::setDT(gwas1); data.table::setDT(gwas2)
  data.table::setDT(LDSC_Tran); data.table::setDT(LDSC1); data.table::setDT(LDSC2)

  # SNP overlap of GWAS only, then use that as the skeleton for all joins
  SNPintersect = intersect(gwas1$SNP, gwas2$SNP)
  skel = data.table::data.table(SNP = SNPintersect)
  data.table::setkey(skel, SNP)

  # Left join in a fixed order, nomatch=0 to keep strict common SNPs
  data.table::setkey(gwas1, SNP);      g1 = gwas1[skel, nomatch = 0]
  data.table::setkey(gwas2, SNP);      g2 = gwas2[skel, nomatch = 0]
  data.table::setkey(LDSC_Tran, SNP);  LDSC_Tran = LDSC_Tran[skel, nomatch = 0]
  data.table::setkey(LDSC1, SNP);      LDSC1  = LDSC1[skel, nomatch = 0]
  data.table::setkey(LDSC2, SNP);      LDSC2  = LDSC2[skel, nomatch = 0]

  # Quick sanity: all same length and same SNP order
  stopifnot(nrow(g1) == nrow(g2), nrow(g1) == nrow(LDSC_Tran),
            nrow(g1) == nrow(LDSC1), nrow(g1) == nrow(LDSC2))
  stopifnot(identical(g1$SNP, g2$SNP),
            identical(g1$SNP, LDSC_Tran$SNP),
            identical(g1$SNP, LDSC1$SNP),
            identical(g1$SNP, LDSC2$SNP))

  M = nrow(g1)

  #----------------------- Construct regressor and weights --------------------------
  # Regression regressor uses cross-ancestry LD score:
  #   LDSC_Tran = LDSC_Tran * sqrt(N1*N2) / M
  z  = g1$Zscore * g2$Zscore
  l  = LDSC_Tran$LDSC * sqrt(g1$N) * sqrt(g2$N) / M

  # Initial WLS with standard LDSC-style weights.
  # Typical bivariate LDSC variance model ≈ (1 + h21 * N1 * LDSC1 / M) * (1 + h22 * N2 * LDSC2 / M)
  # plus a regression mean term (ecov + gcov*LDSC_Tran) that contributes 2 * mean^2 in the denom.
  # Start with product form (recommended); you can switch back to your ratio if you insist.
  X  = cbind(1, l)

  w0_core = (1 + LDSC1$LDSC * g1$N * h21 / M) * (1 + LDSC2$LDSC * g2$N * h22 / M)
  w       = 1 / w0_core

  XtX = CppMatrix::matrixMultiply(t(X), X * w)
  Xty = CppMatrix::matrixVectorMultiply(t(X), z * w)
  beta = as.vector(solve(XtX, Xty))
  ecov.ini = beta[1]; gcov.ini = beta[2]

  # Reweight for efficiency with the mean term included
  mu  = ecov.ini + gcov.ini * l
  w   = 1 / (w0_core + 2 * mu^2)

  XtX = CppMatrix::matrixMultiply(t(X), X * w)
  Xty = CppMatrix::matrixVectorMultiply(t(X), z * w)
  beta = as.vector(solve(XtX, Xty))
  ecov = beta[1]; gcov = beta[2]

  # One more iter (optional)
  mu  = ecov + gcov * l
  w   = 1 / (w0_core + 2 * mu^2)
  XtX = CppMatrix::matrixMultiply(t(X), X * w)
  Xty = CppMatrix::matrixVectorMultiply(t(X), z * w)
  beta = as.vector(solve(XtX, Xty))
  ecov = beta[1]; gcov = beta[2]

  #----------------------- Block bootstrap SEs (same blocks across all tables) ------
  if (sampling.time > 0) {
    blocksize = floor(M / nblock)
    idx_mat   = matrix(seq_len(nblock * blocksize), nrow = blocksize, ncol = nblock)
    gap       = M - nblock * blocksize

    gcov_vec = numeric(sampling.time)
    ecov_vec = numeric(sampling.time)

    for (i in seq_len(sampling.time)) {
      cols = sample.int(nblock, nblock, replace = TRUE)
      ind  = as.vector(idx_mat[, cols]) + if (gap > 0) sample.int(gap, 1) else 0

      X1 = cbind(1, l[ind])
      z1 = z[ind]
      w1 = w[ind]

      XtX = CppMatrix::matrixMultiply(t(X1), X1 * w1)
      Xty = CppMatrix::matrixVectorMultiply(t(X1), z1 * w1)
      b   = as.vector(solve(XtX, Xty))

      ecov_vec[i] = b[1]
      gcov_vec[i] = b[2]
    }

    ecov.se = sqrt(mean((ecov_vec - ecov)^2))
    gcov.se = sqrt(mean((gcov_vec - gcov)^2))
  } else {
    ecov.se = NA_real_
    gcov.se = NA_real_
  }

  data.frame(ecov = ecov, ecov.se = ecov.se,
             gcov = gcov, gcov.se = gcov.se,
             M = M)
}
