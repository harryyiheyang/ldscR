#' Univariate LD Score Regression (LDSC) for SNP-Heritability
#'
#' `ldsc.univ` estimates SNP-heritability for a single trait using
#' Linkage Disequilibrium Score Regression (LDSC). The function harmonizes the
#' GWAS summary statistics with an LD-score table, fits a weighted least squares
#' (WLS) regression of \eqn{Z_j^2} on the LD regressor, and performs an optional
#' block bootstrap for standard errors.
#'
#' @param gwas A `data.frame` or `data.table` of GWAS summary statistics for one trait.
#'   Required columns:
#'   \itemize{
#'     \item `SNP` — SNP identifier
#'     \item `Zscore` — Z-score of the SNP–trait association
#'     \item `N` — per-SNP sample size
#'     \item `A1` — effect allele
#'     \item `A2` — reference allele
#'   }
#'   Only `SNP`, `Zscore`, and `N` are used in the computation; alleles are kept
#'   for consistency checks and downstream use.
#' @param LDSC A `data.frame` or `data.table` containing LD scores.
#'   Required columns:
#'   \itemize{
#'     \item `SNP` — SNP identifier
#'     \item `LDSC` — LD score
#'   }
#' @param nblock Integer. Number of blocks for block bootstrap.
#' @param sampling.time Integer. Number of bootstrap replicates. If `0`, only
#'   point estimates are returned.
#'
#' @details
#' Let \eqn{M} be the number of SNPs after harmonization. The regressor is
#' \eqn{\ell_j = \mathrm{LDSC}_j \, N_j / M}, and we regress \eqn{Z_j^2} on
#' \eqn{(1, \ell_j)} with WLS. Initial weights use a small ridge factor
#' \eqn{(1 + 0.1 \, \ell_j)^{-2}}, followed by an efficiency reweighting that
#' plugs in the estimated heritability \eqn{\hat h^2}:
#' \eqn{w_j = (1 + \ell_j \, \hat h^2)^{-2}}.
#'
#' @return A `data.frame` with:
#' \itemize{
#'   \item `intercept` — estimated LDSC intercept
#'   \item `intercept.se` — standard error of `intercept` (bootstrap; `NA` if no bootstrap)
#'   \item `h2` — estimated SNP-heritability (slope)
#'   \item `h2.se` — standard error of `h2` (bootstrap; `NA` if no bootstrap)
#'   \item `M` — number of SNPs analyzed
#' }
#'
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom data.table setDT setkey
#' @export
ldsc.univ = function(gwas, LDSC, sampling.time = 0, nblock = 500) {
############################# PART 1 #########################################
gwas <- data.table::copy(gwas)
LDSC <- data.table::copy(LDSC)
t0 <- Sys.time()
data.table::setDT(LDSC); data.table::setkey(LDSC, SNP)
skel <- LDSC[, .(SNP)]
data.table::setkey(skel, SNP)
data.table::setDT(gwas); data.table::setkey(gwas, SNP)
gwas <- gwas[skel, nomatch = 0]
gwas <- LDSC[gwas, nomatch = 0]
M <- nrow(gwas)
t0 <- Sys.time() - t0
print("Processing data"); print(t0)

.wls_fit <- function(l, z, w) {
X   <- cbind(1, l)
XtX <- CppMatrix::matrixMultiply(t(X), X * w)
Xty <- CppMatrix::matrixVectorMultiply(t(X), z * w)
as.vector(solve(XtX) %*% Xty)
}

t1 <- Sys.time()
z <- gwas$Zscore^2
l <- gwas$LDSC * gwas$N / M
w <- 1 / pmax(1, gwas$LDSC) / (1 + 0.1 * l)^2
beta <- .wls_fit(l, z, w)
h2.ini <- beta[2]
intercept.ini <- beta[1]
t1 <- Sys.time() - t1
print("Initial Heritability Estimate"); print(t1)
############################# PART 2 #########################################
t2 <- Sys.time()
w  <- 1 / pmax(1, gwas$LDSC) / (1 + l * h2.ini)^2
beta <- .wls_fit(l, z, w)
h2 <- beta[2]; intercept <- beta[1]

w  <- 1 / pmax(1, gwas$LDSC) / (1 + l * h2)^2
beta <- .wls_fit(l, z, w)
h2 <- beta[2]; intercept <- beta[1]
t2 <- Sys.time() - t2
print("Heritability Estimation (Reweighted)"); print(t2)
############################# PART 3 #########################################
if (sampling.time > 0) {
t3 <- Sys.time()

X <- cbind(1, l)
blocksize <- floor(M / nblock)
rem <- M - blocksize * nblock
if (rem >= 0.5 * blocksize) {
block_lengths <- c(rep(blocksize, nblock), rem)
} else {
if (rem > 0) {
block_lengths <- c(rep(blocksize, nblock - 1), blocksize + rem)
} else {
block_lengths <- rep(blocksize, nblock)
}
}
B <- length(block_lengths)
starts <- cumsum(c(1, head(block_lengths, -1)))
IndexList <- lapply(seq_len(B), function(b) starts[b]:(starts[b] + block_lengths[b] - 1))

XtXlist <- vector("list", B)
Xtylist <- vector("list", B)
for (b in seq_len(B)) {
indi <- IndexList[[b]]
Xi   <- X[indi, , drop = FALSE]
wi   <- w[indi]
zi   <- z[indi]
XtXlist[[b]] <- CppMatrix::matrixMultiply(t(Xi), Xi * wi)
Xtylist[[b]] <- CppMatrix::matrixVectorMultiply(t(Xi), zi * wi)
}

h2.vec <- numeric(sampling.time)
intercept.vec <- numeric(sampling.time)

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
beta <- as.vector(solve(XtXi, Xtyi))
intercept.vec[r] <- beta[1]
h2.vec[r]        <- beta[2]
}

t3 <- Sys.time() - t3
print("Standard Error Estimation"); print(t3)

intercept.se <- sd(intercept.vec)
h2.se        <- sd(h2.vec)
} else {
h2.se <- intercept.se <- NA_real_
}

data.frame(intercept = intercept, intercept.se = intercept.se,
h2 = h2, h2.se = h2.se,
M = M)
}
