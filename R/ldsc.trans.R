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
ldsc.trans = function(gwas1, gwas2, h21, h22, LDSC1, LDSC2, LDSC_Tran,sampling.time = 0,nblock = 500) {
gwas1 <- data.table::copy(gwas1)
LDSC1 <- data.table::copy(LDSC1)
gwas2 <- data.table::copy(gwas2)
LDSC2 <- data.table::copy(LDSC2)

data.table::setDT(gwas1); data.table::setDT(gwas2)
data.table::setDT(LDSC1); data.table::setDT(LDSC2); data.table::setDT(LDSC_Tran)

SNPintersect <- Reduce(intersect, list(gwas1$SNP, gwas2$SNP, LDSC_Tran$SNP, LDSC1$SNP, LDSC2$SNP))
skel <- data.table::data.table(SNP = SNPintersect); data.table::setkey(skel, SNP)

data.table::setkey(gwas1, SNP);     g1  <- gwas1[skel, nomatch = 0]
data.table::setkey(gwas2, SNP);     g2  <- gwas2[skel, nomatch = 0]
data.table::setkey(LDSC_Tran, SNP); L12 <- LDSC_Tran[skel, nomatch = 0]
data.table::setkey(LDSC1, SNP);     L1  <- LDSC1[skel, nomatch = 0]
data.table::setkey(LDSC2, SNP);     L2  <- LDSC2[skel, nomatch = 0]

stopifnot(identical(g1$SNP, g2$SNP),
identical(g1$SNP, L12$SNP),
identical(g1$SNP, L1$SNP),
identical(g1$SNP, L2$SNP))

M <- nrow(g1)

z <- g1$Zscore * g2$Zscore
l <- L12$LDSC * sqrt(g1$N) * sqrt(g2$N) / M
X <- cbind(1, l)

base_var <- (1 + L1$LDSC * g1$N * h21 / M) * (1 + L2$LDSC * g2$N * h22 / M)
w <- 1 / (base_var * sqrt(pmax(1, L1$LDSC) * pmax(1, L2$LDSC)))

.fit <- function(X, z, w) {
XtX <- CppMatrix::matrixMultiply(t(X), X * w)
Xty <- CppMatrix::matrixVectorMultiply(t(X), z * w)
as.vector(solve(XtX, Xty))
}

b <- .fit(X, z, w); ecov <- b[1]; gcov <- b[2]
mu <- ecov + gcov * l
w  <- 1 / ((base_var + 2 * mu^2) * sqrt(pmax(1, L1$LDSC) * pmax(1, L2$LDSC)))
b  <- .fit(X, z, w); ecov <- b[1]; gcov <- b[2]

mu <- ecov + gcov * l
w  <- 1 / ((base_var + 2 * mu^2) * sqrt(pmax(1, L1$LDSC) * pmax(1, L2$LDSC)))
b  <- .fit(X, z, w); ecov <- b[1]; gcov <- b[2]

if (sampling.time > 0) {
blocksize <- floor(M / nblock)
rem <- M - blocksize * nblock
if (rem >= 0.5 * blocksize) {
blen <- c(rep(blocksize, nblock), rem)
} else {
blen <- if (rem > 0) c(rep(blocksize, nblock - 1), blocksize + rem) else rep(blocksize, nblock)
}
B <- length(blen)
starts <- cumsum(c(1, head(blen, -1)))
IndexList <- lapply(seq_len(B), function(bi) starts[bi]:(starts[bi] + blen[bi] - 1))

XtXlist <- vector("list", B); Xtylist <- vector("list", B)
for (bi in seq_len(B)) {
idx <- IndexList[[bi]]
Xi  <- X[idx, , drop = FALSE]
wi  <- w[idx]
zi  <- z[idx]
XtXlist[[bi]] <- CppMatrix::matrixMultiply(t(Xi), Xi * wi)
Xtylist[[bi]] <- CppMatrix::matrixVectorMultiply(t(Xi), zi * wi)
}

ecov.vec <- gcov.vec <- numeric(sampling.time)
for (r in seq_len(sampling.time)) {
draw <- sample.int(B, B, replace = TRUE)
cnt  <- tabulate(draw, nbins = B)
sel  <- which(cnt > 0L)
XtXi <- matrix(0, 2, 2); Xtyi <- numeric(2)
for (j in sel) {
cj   <- cnt[j]
XtXi <- XtXi + cj * XtXlist[[j]]
Xtyi <- Xtyi + cj * Xtylist[[j]]
}
br <- as.vector(solve(XtXi, Xtyi))
ecov.vec[r] <- br[1]; gcov.vec[r] <- br[2]
}
ecov.se <- sd(ecov.vec)
gcov.se <- sd(gcov.vec)
} else {
ecov.se <- gcov.se <- NA_real_
}

data.frame(ecov = ecov, ecov.se = ecov.se,
gcov = gcov, gcov.se = gcov.se,
M = M)
}
