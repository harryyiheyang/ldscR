# Shared helpers for the LDSC / SumHer regressions.

# Default truncation of |Z| (ldsc removes chi^2 > max(80, 0.001 N); here the same cap is applied by truncation)
.zmax_default <- function(N) sqrt(max(80, 0.001 * max(N, na.rm = TRUE)))

.truncate_z <- function(z, z.max) sign(z) * pmin(abs(z), z.max)

# Contiguous block index list (same rule as the original ldscR code)
.block_index <- function(M, nblock) {
  nblock <- min(nblock, M)
  blocksize <- floor(M / nblock)
  rem <- M - blocksize * nblock
  if (rem >= 0.5 * blocksize) {
    block_lengths <- c(rep(blocksize, nblock), rem)
  } else if (rem > 0) {
    block_lengths <- c(rep(blocksize, nblock - 1), blocksize + rem)
  } else {
    block_lengths <- rep(blocksize, nblock)
  }
  block_lengths <- block_lengths[block_lengths > 0]
  starts <- cumsum(c(1, utils::head(block_lengths, -1)))
  lapply(seq_along(block_lengths), function(b) starts[b]:(starts[b] + block_lengths[b] - 1))
}

# Bootstrap block counts (one replicate per row). A supplied matrix is reused, so several fits can share draws.
.boot_counts <- function(B, sampling.time, block.counts = NULL) {
  if (!is.null(block.counts)) {
    if (ncol(block.counts) != B) stop("block.counts has ", ncol(block.counts), " columns but there are ", B, " blocks")
    return(block.counts)
  }
  cnt <- matrix(0L, nrow = sampling.time, ncol = B)
  for (r in seq_len(sampling.time)) cnt[r, ] <- tabulate(sample.int(B, B, replace = TRUE), nbins = B)
  colnames(cnt) <- paste0("block", seq_len(B))
  cnt
}

# Per-block weighted cross-products: XtX and Xty for rows idx
.block_xprod <- function(X, y, w, IndexList) {
  lapply(IndexList, function(idx) {
    Xi <- X[idx, , drop = FALSE]
    list(XtX = crossprod(Xi, Xi * w[idx]), Xty = as.vector(crossprod(Xi, y[idx] * w[idx])))
  })
}

.sum_xprod <- function(blocks, cnt) {
  sel <- which(cnt > 0)
  XtX <- Reduce(`+`, Map(function(b, c) c * b$XtX, blocks[sel], cnt[sel]))
  Xty <- Reduce(`+`, Map(function(b, c) c * b$Xty, blocks[sel], cnt[sel]))
  list(XtX = XtX, Xty = Xty)
}

# Two-parameter LDSC fit y = a + b l. With two-step, a comes from the step-1 rows (S1) and b from all rows with a fixed.
# s1: XtX/Xty over S1 rows (2x2, 2); s2: XtX/Xty over all rows.
.fit_ab <- function(s1, s2, two.step) {
  if (!two.step) return(as.vector(solve(s2$XtX, s2$Xty)))
  a <- solve(s1$XtX, s1$Xty)[1]
  b <- (s2$Xty[2] - a * s2$XtX[1, 2]) / s2$XtX[2, 2]
  c(a, b)
}
