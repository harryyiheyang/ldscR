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
#'   estimated (with all \eqn{\tau_c}) from the fit with \eqn{\chi^2} truncated at
#'   `two.step` (all SNPs kept), then the \eqn{\tau_c} from the full statistics with
#'   \eqn{a} fixed. `NULL` (default) = one step.
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
#'   `Partition` (one row per column; every SE from the same replicates):
#'   additive decomposition `h2` = \eqn{	au_c \sum_j a_{jc}} (the columns add up to the total
#'   \eqn{h^2}; `Base` = \eqn{	au_{Base} M}), `Prop_h2` = `h2` / total, `Prop_SNPs` =
#'   \eqn{\sum_j a_{jc} / M}, `tau`, `tau.se`, `tau_z`, `tau_p`; then, for 0/1 columns only, the
#'   ldsc `--overlap-annot` quantities (not additive across overlapping columns):
#'   `h2_overlap` = \eqn{(A^	op A 	au)_c} (h2 of the SNPs in category \eqn{c}),
#'   `Prop_h2_overlap`, `Enrichment_overlap` = `Prop_h2_overlap / Prop_SNPs` and
#'   `Enrichment_overlap_p` (two-sided p of \eqn{h^2_c/M_c - (h^2-h^2_c)/(M-M_c) = 0}),
#'   `replicates` (jackknife / bootstrap coefficient estimates), `SE_method`, and with
#'   `ridge = TRUE` the chosen `lambda` and the cross-validation table `cv`.
#'
#' @export
ldsc.strat = function(gwas, SLDSC, annot, sampling.time = 0, nblock = 500, intercetp = T, intercetp.lower = 0.9, intercetp.upper = 2, two.step = NULL, z.max = NULL, ridge = FALSE, lambda = NULL, lambda.grid = c(0, 10^seq(-6, 2, length.out = 49))) {
  gwas <- as.data.frame(gwas); SLDSC <- as.data.frame(SLDSC); annot <- as.data.frame(annot)
  annot <- annot[!duplicated(annot$SNP), ]
  annot_names <- colnames(annot)[-1]
  Apanel <- as.matrix(annot[, -1, drop = FALSE]); storage.mode(Apanel) <- "double"
  is_bin <- apply(Apanel, 2, function(v) all(v %in% c(0, 1)))

  common <- SLDSC$SNP[SLDSC$SNP %in% gwas$SNP & SLDSC$SNP %in% annot$SNP]
  SLDSC <- SLDSC[match(common, SLDSC$SNP), ]
  gwas <- gwas[match(common, gwas$SNP), ]
  block <- NULL
  if ("Block" %in% colnames(SLDSC)) { block <- SLDSC$Block; SLDSC$Block <- NULL }
  wld <- NULL
  if ("wLD" %in% colnames(SLDSC)) { wld <- SLDSC$wLD; SLDSC$wLD <- NULL }
  L <- as.matrix(SLDSC[, -1, drop = FALSE]); storage.mode(L) <- "double"
  if (ncol(L) != length(annot_names)) stop("SLDSC and annot must have the same LD-score / annotation columns")
  if (is.null(wld)) wld <- L[, 1]
  if (is.null(z.max)) z.max <- .zmax_default(gwas$N)

  n <- nrow(gwas)
  jack <- !is.null(block)
  if (jack) {
    grp <- as.integer(factor(block, levels = unique(block))) - 1L
    cnt <- matrix(0, 0, 0)
  } else {
    B <- .n_groups(n, nblock)
    size <- floor(n / min(nblock, n)); lens <- rep(size, min(nblock, n))
    if (B > length(lens)) lens <- c(lens, n - sum(lens)) else lens[length(lens)] <- lens[length(lens)] + n - sum(lens)
    grp <- rep.int(seq_along(lens) - 1L, lens)
    cnt <- .boot_counts(B, sampling.time)
  }
  lam <- if (is.null(lambda)) NaN else lambda
  r <- ldsc_strat_cpp(as.numeric(gwas$Zscore), as.numeric(gwas$N), L, as.numeric(wld), Apanel, grp, jack, cnt,
                      isTRUE(intercetp), intercetp.lower, intercetp.upper, .two_step_arg(two.step), z.max, 3L,
                      isTRUE(ridge), lam, as.numeric(lambda.grid), 1e-12)
  M <- r$M; acov <- as.vector(r$acov); overlap <- r$overlap
  intercept <- r$b[1]; tau <- r$b[-1]; h2 <- r$h2
  est <- r$est
  if (nrow(est) > 1) {
    B <- nrow(est)
    .se <- if (jack) function(v) sqrt((B - 1) / B * sum((v - mean(v))^2)) else stats::sd
    tau.rep <- est[, -1, drop = FALSE]
    h2.rep <- as.vector(tau.rep %*% acov)
    intercept.se <- if (intercetp) .se(est[, 1]) else 0
    tau.se <- apply(tau.rep, 2, .se); h2.se <- .se(h2.rep)
    prop.se <- apply((tau.rep %*% overlap) / h2.rep, 2, .se)
    se.method <- if (jack) paste0("LOBO jackknife, ", B, " blocks") else paste0("block bootstrap, ", ncol(cnt), " contiguous blocks, ", B, " replicates")
  } else {
    intercept.se <- if (intercetp) NA_real_ else 0
    tau.se <- prop.se <- rep(NA_real_, length(tau)); h2.se <- NA_real_; se.method <- "none"
  }
  Partition <- .partition_table(annot_names, tau, overlap, acov, acov, M, is_bin,
                                reps = if (nrow(est) > 1) est[, -1, drop = FALSE] else NULL, se_fun = if (nrow(est) > 1) .se else NULL)
  list(H2_summary = data.frame(Category = c("Total", "Residual", annot_names[-1]), Estimate = c(h2, tau * acov), SE = c(h2.se, tau.se * acov)),
       Parameter_summary = data.frame(Parameter = c("Intercept", "Residual", annot_names[-1]), Estimate = c(intercept, tau), SE = c(intercept.se, tau.se)),
       Partition = Partition, replicates = if (nrow(est) > 1) est else NULL, SE_method = se.method, n_truncated = r$n_truncated, z.max = z.max,
       M = M, lambda = r$lambda, cv = if (length(r$cv)) data.frame(lambda = lambda.grid, cv = r$cv) else NULL,
       annot_var = as.vector(r$avar))
}
