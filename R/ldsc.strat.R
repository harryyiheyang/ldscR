#' Stratified LD Score Regression on the squared-effect scale
#'
#' `ldsc.strat` fits stratified LDSC with the squared standardized effect
#' \eqn{\hat\beta_j^2 = Z_j^2 / N_j} as the outcome:
#' \deqn{E[\hat\beta_j^2] = \sum_c \tau_c \ell_{jc} + a / N_j .}
#' The LD scores enter without the sample size, and the usual LDSC intercept is
#' the coefficient \eqn{a} of the covariate \eqn{1/N_j} (\eqn{a = 1} without
#' confounding). With the weights below this is the same regression as
#' \eqn{Z_j^2} on \eqn{(1, N_j \ell_{jc})}.
#'
#' @param gwas `data.frame` with columns `SNP`, `Zscore`, `N`.
#' @param SLDSC `data.frame` with column `SNP`, then the LD-score columns; the
#'   first LD-score column is the overall LD score (used in the weights). An
#'   optional column `Block` gives the LD-block index of every SNP.
#' @param annot `data.frame` with column `SNP`, then the annotation columns in
#'   the same order as the LD-score columns (the first is the all-ones column).
#' @param sampling.time Number of block-bootstrap replicates; used only when
#'   `SLDSC` has no `Block` column. `0` returns point estimates only.
#' @param nblock Number of contiguous equal-size SNP blocks for the bootstrap;
#'   used only when `SLDSC` has no `Block` column.
#' @param intercetp If `TRUE` the coefficient \eqn{a} of \eqn{1/N_j} is
#'   estimated within `[intercetp.lower, intercetp.upper]`; if `FALSE` it is
#'   fixed at 1.
#' @param intercetp.lower,intercetp.upper Bounds for \eqn{a}.
#' @param z.max Z-scores are truncated to `[-z.max, z.max]` (not removed).
#' @param ridge If `TRUE`, the \eqn{\tau_c} are estimated by ridge regression
#'   with penalty \eqn{\lambda \sum_c v_c \tau_c^2}, where \eqn{v_c} is the
#'   variance of annotation column \eqn{c} over the SNPs (\eqn{p_c(1-p_c)} for a
#'   0/1 column). This penalises the standardized effects
#'   \eqn{\tau_c \, \mathrm{sd}(a_c)} equally. The all-ones column
#'   (\eqn{v = 0}) and the \eqn{1/N} term are not penalised.
#' @param lambda Ridge parameter on the scale of `lambda.grid` (a multiple of
#'   the median of \eqn{(X^\top W X)_{cc} / v_c}). `NULL` chooses it by
#'   leave-one-block-out cross-validation.
#' @param lambda.grid Candidate values for the cross-validation.
#'
#' @details
#' Standard errors. If `SLDSC` has a `Block` column, they are leave-one-block-out
#' (LOBO) jackknife standard errors over the LD blocks, and `sampling.time` and
#' `nblock` are ignored. Otherwise they are block-bootstrap standard errors over
#' `nblock` contiguous blocks (rows must be in genomic order). The bound on
#' \eqn{a} is applied in every leave-out or bootstrap fit, as in the full fit.
#'
#' Weights: \eqn{w_j = 1 / \max(1, \ell_{j1}) / (1/N_j + \hat h^2 \ell_{j1} / M)^2},
#' three reweighting rounds, then fixed.
#'
#' @return A list with `H2_summary` (total and per-column heritability
#'   \eqn{\tau_c \sum_j a_{jc}}), `Parameter_summary` (row `Intercept` = \eqn{a},
#'   then \eqn{\tau_c}), `SE_method`, and with `ridge = TRUE` the chosen `lambda`
#'   and the cross-validation table `cv` (error minus its minimum). Standard
#'   errors under ridge are for the fixed chosen `lambda`.
#'
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @export
ldsc.strat = function(gwas, SLDSC, annot, sampling.time = 0, nblock = 500, intercetp = T, intercetp.lower = 0.9, intercetp.upper = 2, z.max = 50, ridge = FALSE, lambda = NULL, lambda.grid = c(0, 10^seq(-6, 2, length.out = 49))) {
  t0 <- Sys.time()
  gwas <- data.table::copy(as.data.frame(gwas))
  SLDSC <- data.table::copy(as.data.frame(SLDSC))
  annot = data.table::copy(as.data.frame(annot))
  data.table::setDT(SLDSC)
  data.table::setDT(gwas)
  data.table::setDT(annot)

  common_snps <- SLDSC$SNP[SLDSC$SNP %in% gwas$SNP & SLDSC$SNP %in% annot$SNP]
  SLDSC <- SLDSC[match(common_snps, SLDSC$SNP), ]
  gwas <- gwas[match(common_snps, gwas$SNP), ]
  annot <- annot[match(common_snps, annot$SNP), ]
  M <- nrow(gwas)

  block <- NULL
  if ("Block" %in% colnames(SLDSC)) {
    block <- SLDSC$Block
    SLDSC[, Block := NULL]
  }
  t0 <- Sys.time() - t0
  print("Processing data"); print(t0)

  X = as.matrix(SLDSC[, -1])
  Annot = as.matrix(annot[, -1])
  annot_names = colnames(annot)[-1]
  Annot_cov = colSums(Annot)
  Annot_var = pmax(colMeans(Annot^2) - colMeans(Annot)^2, 0)
  remove(SLDSC, annot)

  # solve on unit-diagonal scale (the 1/N column and the LD-score columns differ by many orders of magnitude)
  .solve <- function(A, b) {
    d <- 1 / sqrt(diag(A))
    as.vector(d * qr.solve(A * tcrossprod(d), b * d, tol = 1e-12))
  }
  # first row / column = the 1/N term; returns c(a, tau)
  .fit <- function(XtX, Xty) {
    if (!intercetp) return(c(1, .solve(XtX[-1, -1, drop = FALSE], Xty[-1] - XtX[-1, 1])))
    res <- .solve(XtX, Xty)
    if (res[1] < intercetp.lower) {
      res <- c(intercetp.lower, .solve(XtX[-1, -1, drop = FALSE], Xty[-1] - intercetp.lower * XtX[-1, 1]))
    } else if (res[1] > intercetp.upper) {
      res <- c(intercetp.upper, .solve(XtX[-1, -1, drop = FALSE], Xty[-1] - intercetp.upper * XtX[-1, 1]))
    }
    res
  }

  t1 <- Sys.time()
  zt <- pmin(abs(gwas$Zscore), z.max)
  n.trunc <- sum(abs(gwas$Zscore) > z.max)
  invN <- 1 / gwas$N
  y <- zt^2 * invN
  l1 <- X[, 1]
  lw <- 1 / pmax(1, l1)
  XX <- cbind(invN, X)
  remove(X)

  w <- lw / (invN + 0.1 * l1 / M)^2
  beta_ini <- .fit(CppMatrix::matrixMultiply(t(XX), XX * w), CppMatrix::matrixVectorMultiply(t(XX), y * w))
  h2.ini <- sum(beta_ini[-1] * Annot_cov)
  t1 <- Sys.time() - t1
  print("Initial Heritability Estimate"); print(t1)

  t2 <- Sys.time()
  h2_curr <- h2.ini
  for (i in 1:3) {
    w <- lw / (invN + l1 * h2_curr / M)^2
    beta_tmp <- .fit(CppMatrix::matrixMultiply(t(XX), XX * w), CppMatrix::matrixVectorMultiply(t(XX), y * w))
    intercept <- beta_tmp[1]
    tau <- beta_tmp[-1]
    h2_curr <- sum(tau * Annot_cov)
  }
  h2 <- h2_curr
  t2 <- Sys.time() - t2
  print("Heritability Estimation (Reweighted)"); print(t2)

  K <- ncol(XX)
  se.method <- "none"
  Pen <- matrix(0, K, K)
  lambda.used <- NA_real_
  cv.table <- NULL
  if (!is.null(block) || sampling.time > 0 || ridge) {
    t3 <- Sys.time()
    if (!is.null(block)) {
      IndexList <- split(seq_len(M), block)
    } else {
      blocksize <- floor(M / nblock)
      rem <- M - blocksize * nblock
      if (rem >= 0.5 * blocksize) {
        block_lengths <- c(rep(blocksize, nblock), rem)
      } else if (rem > 0) {
        block_lengths <- c(rep(blocksize, nblock - 1), blocksize + rem)
      } else {
        block_lengths <- rep(blocksize, nblock)
      }
      starts <- cumsum(c(1, head(block_lengths, -1)))
      IndexList <- lapply(seq_along(block_lengths), function(b) starts[b]:(starts[b] + block_lengths[b] - 1))
    }
    B <- length(IndexList)
    XtXlist <- vector("list", B)
    Xtylist <- vector("list", B)
    for (b in seq_len(B)) {
      idx <- IndexList[[b]]
      Xi <- XX[idx, , drop = FALSE]
      wi <- w[idx]
      XtXlist[[b]] <- CppMatrix::matrixMultiply(t(Xi), Xi * wi)
      Xtylist[[b]] <- CppMatrix::matrixVectorMultiply(t(Xi), y[idx] * wi)
    }
    XtX_all <- Reduce(`+`, XtXlist)
    Xty_all <- Reduce(`+`, Xtylist)

    if (ridge) {
      # penalty lambda * sum_c v_c tau_c^2; the 1/N term (first) and columns with v = 0 are free
      pen.d <- c(0, Annot_var)
      lam.unit <- stats::median((diag(XtX_all) / pen.d)[pen.d > 0])
      if (is.null(lambda)) {
        cv.err <- numeric(length(lambda.grid))
        for (g in seq_along(lambda.grid)) {
          Pg <- diag(lambda.grid[g] * lam.unit * pen.d)
          e <- 0
          for (b in seq_len(B)) {
            bb <- .fit(XtX_all + Pg - XtXlist[[b]], Xty_all - Xtylist[[b]])
            e <- e - 2 * sum(bb * Xtylist[[b]]) + sum(bb * (XtXlist[[b]] %*% bb))
          }
          cv.err[g] <- e
        }
        cv.table <- data.frame(lambda = lambda.grid, cv = cv.err - min(cv.err))
        lambda.used <- lambda.grid[which.min(cv.err)]
      } else {
        lambda.used <- lambda
      }
      Pen <- diag(lambda.used * lam.unit * pen.d)
      beta_tmp <- .fit(XtX_all + Pen, Xty_all)
      intercept <- beta_tmp[1]
      tau <- beta_tmp[-1]
      h2 <- sum(tau * Annot_cov)
    }

    if (!is.null(block)) {
      # leave-one-block-out jackknife
      se.method <- paste0("LOBO jackknife, ", B, " blocks")
      est.mat <- matrix(NA_real_, B, K)
      for (b in seq_len(B)) est.mat[b, ] <- .fit(XtX_all + Pen - XtXlist[[b]], Xty_all - Xtylist[[b]])
      h2.vec <- as.vector(est.mat[, -1, drop = FALSE] %*% Annot_cov)
      .jk <- function(v) sqrt((B - 1) / B * sum((v - mean(v))^2))
      intercept.se <- .jk(est.mat[, 1])
      tau.se <- apply(est.mat[, -1, drop = FALSE], 2, .jk)
      h2.se <- .jk(h2.vec)
    } else if (sampling.time > 0) {
      se.method <- paste0("block bootstrap, ", B, " contiguous blocks, ", sampling.time, " replicates")
      h2.vec <- numeric(sampling.time)
      intercept.vec <- numeric(sampling.time)
      tau.mat <- matrix(NA, nrow = sampling.time, ncol = length(tau))
      for (r in seq_len(sampling.time)) {
        cnt <- tabulate(sample.int(B, B, replace = TRUE), nbins = B)
        XtXi_sum <- matrix(0, K, K)
        Xtyi_sum <- numeric(K)
        for (j in which(cnt > 0L)) {
          XtXi_sum <- XtXi_sum + cnt[j] * XtXlist[[j]]
          Xtyi_sum <- Xtyi_sum + cnt[j] * Xtylist[[j]]
        }
        res_beta <- .fit(XtXi_sum + Pen, Xtyi_sum)
        intercept.vec[r] <- res_beta[1]
        tau.mat[r, ] <- res_beta[-1]
        h2.vec[r] <- sum(res_beta[-1] * Annot_cov)
      }
      intercept.se <- sd(intercept.vec)
      h2.se <- sd(h2.vec)
      tau.se <- apply(tau.mat, 2, sd)
    } else {
      h2.se <- NA_real_
      intercept.se <- NA_real_
      tau.se <- rep(NA_real_, length(tau))
    }
    if (!intercetp) intercept.se <- 0
    t3 <- Sys.time() - t3
    print("Standard Error Estimation"); print(t3)
  } else {
    h2.se <- NA_real_
    intercept.se <- if (intercetp) NA_real_ else 0
    tau.se <- rep(NA_real_, length(tau))
  }

  tau_h2 <- tau * Annot_cov
  tau_h2.se <- tau.se * Annot_cov

  H2_summary <- data.frame(
    Category = c("Total", "Residual", annot_names[-1]),
    Estimate = c(h2, tau_h2),
    SE = c(h2.se, tau_h2.se)
  )

  Parameter_summary <- data.frame(
    Parameter = c("Intercept", "Residual", annot_names[-1]),
    Estimate = c(intercept, tau),
    SE = c(intercept.se, tau.se)
  )

  list(H2_summary = H2_summary, Parameter_summary = Parameter_summary, SE_method = se.method, n_truncated = n.trunc,
       lambda = lambda.used, cv = cv.table, annot_var = Annot_var)
}
