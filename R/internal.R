# Shared helpers (the numerical work is in src/).

# Default truncation of |Z| (ldsc removes chi^2 > max(80, 0.001 N); here the same cap is applied by truncation)
.zmax_default <- function(N) sqrt(max(80, 0.001 * max(N, na.rm = TRUE)))

# Number of contiguous bootstrap groups for M SNPs (same rule as equal_groups in src/ldsc_kernel.h)
.n_groups <- function(M, nblock) {
  nblock <- min(nblock, M)
  size <- floor(M / nblock); rem <- M - size * nblock
  nblock + (rem >= 0.5 * size && rem > 0)
}

# Bootstrap block counts (replicates x groups); drawn in R so that set.seed() controls them
.boot_counts <- function(B, sampling.time) {
  cnt <- matrix(0, nrow = sampling.time, ncol = B)
  for (r in seq_len(sampling.time)) cnt[r, ] <- tabulate(sample.int(B, B, replace = TRUE), nbins = B)
  colnames(cnt) <- paste0("block", seq_len(B))
  cnt
}

.two_step_arg <- function(two.step) if (is.null(two.step)) -1 else two.step

.ldsc_bootstrap_sd <- function(est) if (nrow(est) > 1) apply(est, 2, stats::sd) else rep(NA_real_, ncol(est))
