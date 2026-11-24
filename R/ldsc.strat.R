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
#'   optional column `Block` gives the LD-block index of every SNP, and an
#'   optional column `wLD` the LD score for the over-counting weight
#'   (default: the first LD-score column).
#' @param annot `data.frame` with column `SNP`, then the annotation columns in
#'   the same order as the LD-score columns (the first is the all-ones column).
#'   It must cover the whole panel (every SNP the LD scores were summed over),
#'   because \eqn{h^2_c = \tau_c \sum_{j \in panel} a_{jc}}.
#' @param sampling.time Number of block-bootstrap replicates; used only when
#'   `SLDSC` has no `Block` column. `0` returns point estimates only.
#' @param nblock Number of contiguous equal-size SNP blocks for the bootstrap;
#'   used only when `SLDSC` has no `Block` column.
#' @param intercetp If `TRUE` the coefficient \eqn{a} of \eqn{1/N_j} is
#'   estimated within `[intercetp.lower, intercetp.upper]`; if `FALSE` it is
#'   fixed at 1.
#' @param intercetp.lower,intercetp.upper Bounds for \eqn{a}.
#' @param two.step Chi-square cutoff for the two-step estimator: \eqn{a} is
#'   estimated (with all \eqn{\tau_c}) from SNPs with \eqn{\chi^2 <} `two.step`,
#'   then the \eqn{\tau_c} from all SNPs with \eqn{a} fixed. `NULL` (default) =
#'   one step: selecting SNPs on \eqn{\chi^2} biases the step-1 intercept upward
#'   when many SNPs exceed the cutoff, and ldsc does not allow the two-step with
#'   partitioned LD scores.
#' @param z.max Z-scores are truncated to `[-z.max, z.max]` (not removed).
#'   `NULL` uses \eqn{\sqrt{\max(80, 0.001 N_{max})}}.
#' @param ridge If `TRUE`, the \eqn{\tau_c} are estimated by ridge regression
#'   with penalty \eqn{\lambda \sum_c v_c \tau_c^2}, where \eqn{v_c} is the
#'   variance of annotation column \eqn{c} over the SNPs. The all-ones column
#'   (\eqn{v = 0}) and the \eqn{1/N} term are not penalised.
#' @param lambda Ridge parameter on the scale of `lambda.grid` (a multiple of
#'   the median of \eqn{(X^\top W X)_{cc} / v_c}). `NULL` chooses it by
#'   leave-one-block-out cross-validation.
#' @param lambda.grid Candidate values for the cross-validation.
#'
#' @details
#' Standard errors. If `SLDSC` has a `Block` column, they are leave-one-block-out
#' (LOBO) jackknife standard errors over the LD blocks. Otherwise they are
#' block-bootstrap standard errors over `nblock` contiguous blocks (rows must
#' be in genomic order). The bound on \eqn{a} and the two-step are applied in
#' every leave-out or bootstrap fit, as in the full fit.
#'
#' Weights: \eqn{w_j = 1 / \max(1, \ell^{w}_j) / (1/N_j + \hat h^2 \ell_{j1} / M)^2},
#' with \eqn{M} the number of panel SNPs; three reweighting rounds, then fixed.
#'
#' @return A list with `H2_summary` (total and per-column \eqn{\tau_c \sum_j a_{jc}}),
#'   `Parameter_summary` (row `Intercept` = \eqn{a}, then \eqn{\tau_c}),
#'   `Enrichment_summary` (for 0/1 annotations: share of SNPs, share of
#'   \eqn{h^2} carried by the SNPs in the category including overlaps,
#'   \eqn{(A^\top A \tau)_c / h^2}, and their ratio), `SE_method`, and with
#'   `ridge = TRUE` the chosen `lambda` and the cross-validation table `cv`.
#'
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @export
ldsc.strat = function(gwas, SLDSC, annot, sampling.time = 0, nblock = 500, intercetp = T, intercetp.lower = 0.9, intercetp.upper = 2, two.step = NULL, z.max = NULL, ridge = FALSE, lambda = NULL, lambda.grid = c(0, 10^seq(-6, 2, length.out = 49))) {
  t0 <- Sys.time()
  gwas <- as.data.frame(gwas)
  SLDSC <- as.data.frame(SLDSC)
  annot <- as.data.frame(annot)
  annot <- annot[!duplicated(annot$SNP), ]

  # panel-level quantities (before restricting to the GWAS SNPs)
  Annot_full <- as.matrix(annot[, -1, drop = FALSE])
  annot_names <- colnames(annot)[-1]
  M <- nrow(Annot_full)
  Annot_cov <- colSums(Annot_full)
  Annot_var <- pmax(colMeans(Annot_full^2) - colMeans(Annot_full)^2, 0)
  is_bin <- apply(Annot_full, 2, function(v) all(v %in% c(0, 1)))
  Overlap <- crossprod(Annot_full)
  remove(Annot_full)

  common_snps <- SLDSC$SNP[SLDSC$SNP %in% gwas$SNP & SLDSC$SNP %in% annot$SNP]
  SLDSC <- SLDSC[match(common_snps, SLDSC$SNP), ]
  gwas <- gwas[match(common_snps, gwas$SNP), ]
  M.reg <- nrow(gwas)
  remove(annot)

  block <- NULL
  if ("Block" %in% colnames(SLDSC)) { block <- SLDSC$Block; SLDSC$Block <- NULL }
  wld <- NULL
  if ("wLD" %in% colnames(SLDSC)) { wld <- SLDSC$wLD; SLDSC$wLD <- NULL }
  X <- as.matrix(SLDSC[, -1, drop = FALSE])
  if (ncol(X) != length(annot_names)) stop("SLDSC and annot must have the same LD-score / annotation columns")
  remove(SLDSC)
  t0 <- Sys.time() - t0
  print("Processing data"); print(t0)

  # solve on unit-diagonal scale (the 1/N column and the LD-score columns differ by many orders of magnitude)
  .solve <- function(A, b) {
    d <- 1 / sqrt(diag(A))
    as.vector(d * qr.solve(A * tcrossprod(d), b * d, tol = 1e-12))
  }
  .fix <- function(XtX, Xty, a) c(a, .solve(XtX[-1, -1, drop = FALSE], Xty[-1] - a * XtX[-1, 1]))
  # one-step fit; first row / column = the 1/N term; returns c(a, tau)
  .fit1 <- function(XtX, Xty) {
    if (!intercetp) return(.fix(XtX, Xty, 1))
    res <- .solve(XtX, Xty)
    if (res[1] < intercetp.lower) res <- .fix(XtX, Xty, intercetp.lower)
    if (res[1] > intercetp.upper) res <- .fix(XtX, Xty, intercetp.upper)
    res
  }
  use2 <- intercetp && !is.null(two.step)
  # s1 = sums over step-1 SNPs, s2 = sums over all SNPs
  .fit <- function(XtX1, Xty1, XtX2, Xty2) {
    if (!use2) return(.fit1(XtX2, Xty2))
    .fix(XtX2, Xty2, .fit1(XtX1, Xty1)[1])
  }

  t1 <- Sys.time()
  if (is.null(z.max)) z.max <- .zmax_default(gwas$N)
  n.trunc <- sum(abs(gwas$Zscore) > z.max)
  zt <- pmin(abs(gwas$Zscore), z.max)
  invN <- 1 / gwas$N
  y <- zt^2 * invN
  S1 <- if (use2) zt^2 < two.step else rep(TRUE, M.reg)
  l1 <- X[, 1]
  lw <- 1 / pmax(1, if (is.null(wld)) l1 else wld)
  XX <- cbind(invN, X)
  remove(X)
  .xp <- function(w) {
    w1 <- w * S1
    list(CppMatrix::matrixMultiply(t(XX), XX * w1), CppMatrix::matrixVectorMultiply(t(XX), y * w1),
         CppMatrix::matrixMultiply(t(XX), XX * w), CppMatrix::matrixVectorMultiply(t(XX), y * w))
  }

  w <- lw / (invN + 0.1 * l1 / M)^2
  beta_ini <- do.call(.fit, .xp(w))
  h2_curr <- sum(beta_ini[-1] * Annot_cov)
  t1 <- Sys.time() - t1
  print("Initial Heritability Estimate"); print(t1)

  t2 <- Sys.time()
  for (i in 1:3) {
    w <- lw / (invN + l1 * max(h2_curr, 0) / M)^2
    beta_tmp <- do.call(.fit, .xp(w))
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
  est.mat <- NULL
  if (!is.null(block) || sampling.time > 0 || ridge) {
    t3 <- Sys.time()
    IndexList <- if (!is.null(block)) split(seq_len(M.reg), block) else .block_index(M.reg, nblock)
    B <- length(IndexList)
    bl1 <- .block_xprod(XX, y, w * S1, IndexList)
    bl2 <- .block_xprod(XX, y, w, IndexList)
    all1 <- .sum_xprod(bl1, rep(1L, B)); all2 <- .sum_xprod(bl2, rep(1L, B))
    .leave <- function(b, P) .fit(all1$XtX + P - bl1[[b]]$XtX, all1$Xty - bl1[[b]]$Xty,
                                  all2$XtX + P - bl2[[b]]$XtX, all2$Xty - bl2[[b]]$Xty)

    if (ridge) {
      # penalty lambda * sum_c v_c tau_c^2; the 1/N term (first) and columns with v = 0 are free
      pen.d <- c(0, Annot_var)
      lam.unit <- stats::median((diag(all2$XtX) / pen.d)[pen.d > 0])
      if (is.null(lambda)) {
        cv.err <- numeric(length(lambda.grid))
        for (g in seq_along(lambda.grid)) {
          Pg <- diag(lambda.grid[g] * lam.unit * pen.d)
          e <- 0
          for (b in seq_len(B)) {
            bb <- .leave(b, Pg)
            e <- e - 2 * sum(bb * bl2[[b]]$Xty) + sum(bb * (bl2[[b]]$XtX %*% bb))
          }
          cv.err[g] <- e
        }
        cv.table <- data.frame(lambda = lambda.grid, cv = cv.err - min(cv.err))
        lambda.used <- lambda.grid[which.min(cv.err)]
      } else {
        lambda.used <- lambda
      }
      Pen <- diag(lambda.used * lam.unit * pen.d)
      beta_tmp <- .fit(all1$XtX + Pen, all1$Xty, all2$XtX + Pen, all2$Xty)
      intercept <- beta_tmp[1]
      tau <- beta_tmp[-1]
      h2 <- sum(tau * Annot_cov)
    }

    if (!is.null(block)) {
      se.method <- paste0("LOBO jackknife, ", B, " blocks")
      est.mat <- t(sapply(seq_len(B), function(b) .leave(b, Pen)))
      .se <- function(v) sqrt((B - 1) / B * sum((v - mean(v))^2))
    } else if (sampling.time > 0) {
      se.method <- paste0("block bootstrap, ", B, " contiguous blocks, ", sampling.time, " replicates")
      cnt <- .boot_counts(B, sampling.time)
      est.mat <- t(apply(cnt, 1, function(cc) {
        s1 <- .sum_xprod(bl1, cc); s2 <- .sum_xprod(bl2, cc)
        .fit(s1$XtX + Pen, s1$Xty, s2$XtX + Pen, s2$Xty)
      }))
      .se <- function(v) sd(v)
    }
    t3 <- Sys.time() - t3
    print("Standard Error Estimation"); print(t3)
  }

  tau_h2 <- tau * Annot_cov
  h2cat <- as.vector(Overlap %*% tau)
  if (!is.null(est.mat)) {
    tau.rep <- est.mat[, -1, drop = FALSE]
    h2.rep <- as.vector(tau.rep %*% Annot_cov)
    intercept.se <- if (intercetp) .se(est.mat[, 1]) else 0
    tau.se <- apply(tau.rep, 2, .se)
    h2.se <- .se(h2.rep)
    prop.rep <- (tau.rep %*% Overlap) / h2.rep
    prop.se <- apply(prop.rep, 2, .se)
  } else {
    intercept.se <- if (intercetp) NA_real_ else 0
    tau.se <- rep(NA_real_, length(tau)); h2.se <- NA_real_; prop.se <- rep(NA_real_, length(tau))
  }

  H2_summary <- data.frame(
    Category = c("Total", "Residual", annot_names[-1]),
    Estimate = c(h2, tau_h2),
    SE = c(h2.se, tau.se * Annot_cov)
  )
  Parameter_summary <- data.frame(
    Parameter = c("Intercept", "Residual", annot_names[-1]),
    Estimate = c(intercept, tau),
    SE = c(intercept.se, tau.se)
  )
  prop_snps <- Annot_cov / M
  prop_h2 <- h2cat / h2
  Enrichment_summary <- data.frame(
    Category = annot_names, Prop_SNPs = prop_snps, Prop_h2 = prop_h2, Prop_h2.se = prop.se,
    Enrichment = prop_h2 / prop_snps, Enrichment.se = prop.se / prop_snps
  )[is_bin, , drop = FALSE]
  rownames(Enrichment_summary) <- NULL

  list(H2_summary = H2_summary, Parameter_summary = Parameter_summary, Enrichment_summary = Enrichment_summary,
       SE_method = se.method, n_truncated = n.trunc, z.max = z.max, M = M, M.reg = M.reg,
       lambda = lambda.used, cv = cv.table, annot_var = Annot_var)
}
