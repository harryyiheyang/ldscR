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
ldsc.bicov = function(gwas1, gwas2, h21, h22, LDSC, sampling.time = 0, nblock = 500) {
  gwas1 <- data.table::copy(gwas1)
  LDSC <- data.table::copy(LDSC)
  gwas2 <- data.table::copy(gwas2)

data.table::setDT(LDSC); data.table::setkey(LDSC, SNP)
data.table::setDT(gwas1); data.table::setDT(gwas2)

t0 <- Sys.time()
SNPintersect <- intersect(intersect(gwas1$SNP, gwas2$SNP), LDSC$SNP)
skel <- data.table::data.table(SNP = SNPintersect); data.table::setkey(skel, SNP)

#
gwas1 <- gwas1[skel, on = "SNP", nomatch = 0]
gwas2 <- gwas2[skel, on = "SNP", nomatch = 0]
gwas1 <- LDSC[gwas1, on = "SNP", nomatch = 0]
gwas2 <- LDSC[gwas2, on = "SNP", nomatch = 0]

stopifnot(identical(gwas1$SNP, gwas2$SNP))
M <- nrow(gwas1)
t0 <- Sys.time() - t0
print("Processing data"); print(t0)

.wls_fit <- function(l, z, w) {
X <- cbind(1, l)
XtX <- CppMatrix::matrixMultiply(t(X), X * w)
Xty <- CppMatrix::matrixVectorMultiply(t(X), z * w)
as.vector(solve(XtX, Xty))
}

t1 <- Sys.time()
z <- gwas1$Zscore * gwas2$Zscore
l <- gwas1$LDSC * sqrt(gwas1$N) * sqrt(gwas2$N) / M
base_var <- (1 + gwas1$LDSC * gwas1$N / M * h21) * (1 + gwas2$LDSC * gwas2$N / M * h22)
w <- 1 / pmax(1, gwas1$LDSC) / base_var
b <- .wls_fit(l, z, w)
gcov <- b[2]; ecov <- b[1]
t1 <- Sys.time() - t1
print("Initial Genetic Covariance Estimate"); print(t1)

t2 <- Sys.time()
mu <- ecov + gcov * l
w  <- 1 / pmax(1, gwas1$LDSC) / (base_var + 2 * mu^2)
b <- .wls_fit(l, z, w)
gcov <- b[2]; ecov <- b[1]

mu <- ecov + gcov * l
w  <- 1 / pmax(1, gwas1$LDSC) / (base_var + 2 * mu^2)
b <- .wls_fit(l, z, w)
gcov <- b[2]; ecov <- b[1]
t2 <- Sys.time() - t2
print("Genetic Covariance Reweighting"); print(t2)

rg <- gcov / sqrt(h21 * h22)

if (sampling.time > 0) {
t3 <- Sys.time()
X <- cbind(1, l)

blocksize <- floor(M / nblock)
rem <- M - blocksize * nblock
if (rem >= 0.5 * blocksize) {
block_lengths <- c(rep(blocksize, nblock), rem)
} else {
block_lengths <- if (rem > 0) c(rep(blocksize, nblock - 1), blocksize + rem) else rep(blocksize, nblock)
}
B <- length(block_lengths)
starts <- cumsum(c(1, head(block_lengths, -1)))
IndexList <- lapply(seq_len(B), function(bi) starts[bi]:(starts[bi] + block_lengths[bi] - 1))

XtXlist <- vector("list", B)
Xtylist <- vector("list", B)
for (bi in seq_len(B)) {
idx <- IndexList[[bi]]
Xi  <- X[idx, , drop = FALSE]
wi  <- w[idx]
zi  <- z[idx]
XtXlist[[bi]] <- CppMatrix::matrixMultiply(t(Xi), Xi * wi)
Xtylist[[bi]] <- CppMatrix::matrixVectorMultiply(t(Xi), zi * wi)
}

gcov.vec <- ecov.vec <- rg.vec <- numeric(sampling.time)
for (r in seq_len(sampling.time)) {
draw <- sample.int(B, B, replace = TRUE)
cnt  <- tabulate(draw, nbins = B)
sel  <- which(cnt > 0L)

XtXi <- matrix(0, 2, 2)
Xtyi <- numeric(2)
for (j in sel) {
cj   <- cnt[j]
XtXi <- XtXi + cj * XtXlist[[j]]
Xtyi <- Xtyi + cj * Xtylist[[j]]
}
br <- as.vector(solve(XtXi, Xtyi))
ecov.vec[r] <- br[1]
gcov.vec[r] <- br[2]
}

t3 <- Sys.time() - t3
print("Standard Error Estimation"); print(t3)

ecov.se <-sd(ecov.vec)
gcov.se <-sd(gcov.vec)
} else {
ecov.se <- gcov.se <- NA_real_
}

data.frame(
ecov = ecov, ecov.se = ecov.se,
gcov = gcov, gcov.se = gcov.se,
M = M
)
}
