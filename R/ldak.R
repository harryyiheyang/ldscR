#' Read one SBayesRC-format LD eigen block
#'
#' Layout: int32 `m`, int32 `k`, float32 `sumPosEig`, float32 `cutoff`,
#' float32 `lambda[k]`, float32 `U[k][m]` (each eigenvector stored
#' contiguously), so that \eqn{R \approx U \mathrm{diag}(\lambda) U^\top}.
#'
#' @param file Path to `blockX.eigen.bin`.
#' @return A list with `m`, `k`, `sumPosEig`, `cutoff`, `lambda` and `U` (m x k).
#' @export
read_eigen_bin <- function(file) {
  con <- file(file, "rb"); on.exit(close(con))
  mk <- readBin(con, "integer", 2, size = 4)
  hdr <- readBin(con, "numeric", 2, size = 4)
  lam <- readBin(con, "numeric", mk[2], size = 4)
  U <- matrix(readBin(con, "numeric", mk[1] * mk[2], size = 4), nrow = mk[1], ncol = mk[2])
  list(m = mk[1], k = mk[2], sumPosEig = hdr[1], cutoff = hdr[2], lambda = lam, U = U)
}

#' Write one SBayesRC-format LD eigen block
#'
#' Inverse of [read_eigen_bin()]; mainly for tests and for turning a full LD
#' matrix into the eigen format.
#' @param file Output path.
#' @param U Eigenvectors (m x k).
#' @param lambda Eigenvalues (length k).
#' @param cutoff Proportion of variance retained (stored in the header only).
#' @export
write_eigen_bin <- function(file, U, lambda, cutoff = 1) {
  con <- file(file, "wb"); on.exit(close(con))
  writeBin(as.integer(c(nrow(U), ncol(U))), con, size = 4)
  writeBin(c(sum(lambda), cutoff), con, size = 4)
  writeBin(as.numeric(lambda), con, size = 4)
  writeBin(as.numeric(U), con, size = 4)
  invisible(file)
}

# r^2 between rows idx and all SNPs of one block, from X = U sqrt(lambda), rescaled to unit diagonal
.block_r2 <- function(X, d, idx) {
  R <- tcrossprod(X[idx, , drop = FALSE], X)
  R <- R / sqrt(d[idx]) / rep(sqrt(d), each = length(idx))
  R^2
}

#' LDAK heritability-model taggings from an eigen LD panel
#'
#' Computes, for every regression SNP \eqn{j}, the quantities SumHer regresses
#' on, using the within-block LD of an SBayesRC-format eigen panel:
#' \deqn{v_{jc} = \sum_{k} r^2_{jk}\, a_{kc}\, w_k\, [2 f_k (1 - f_k)]^{1 + \alpha}, \qquad
#' t_j = \sum_{k \in reg} r^2_{jk},}
#' where \eqn{w_k} are LDAK weights, \eqn{f_k} the allele frequency,
#' \eqn{a_{kc}} the annotations, and off-diagonal \eqn{r^2} are bias-corrected
#' as \eqn{r^2 - (1 - r^2)/(n_{ref} - 2)} (diagonal = 1). These are the
#' `--calc-tagging` quantities of LDAK, with LD restricted to the panel blocks
#' instead of a 1 Mb window.
#'
#' Weights: `"thin"` is the LDAK-Thin model (`--thin --window-prune`): within
#' each block, SNPs are scanned in order and a SNP is dropped (weight 0) when a
#' later, still-kept SNP has \eqn{r^2 \ge} `thin.r2`; kept SNPs get weight 1.
#' `"none"` sets all weights to 1 (with \eqn{\alpha = -1}: the LDSC/GCTA model).
#' A numeric vector gives the weights directly.
#'
#' @param eigen.files Paths of the eigen blocks, in panel order.
#' @param snp.info `data.frame` with columns `SNP`, `A1`, `A2`, `Freq` (A1
#'   frequency), one row per panel SNP, in the order of the blocks.
#' @param annot Optional `data.frame` with `SNP` then annotation columns (any
#'   SNP not listed gets 0). A `Base` column of ones is always added first.
#' @param alpha Power \eqn{\alpha} (LDAK `--power`): -0.25 for the LDAK-Thin /
#'   Human Default model, -1 for GCTA/LDSC.
#' @param weights `"thin"`, `"none"` or a numeric vector over panel SNPs.
#' @param thin.r2 Pruning threshold for `"thin"`.
#' @param n.ref Sample size of the LD panel for the \eqn{r^2} bias correction
#'   (`NULL` = no correction).
#' @param regression.snps SNPs used as regression SNPs (`NULL` = all panel SNPs).
#'   Restrict this to the SNPs the GWAS will provide.
#' @param chunk Rows of the LD matrix formed at a time.
#' @return A list of class `ldak_tagging`: `tagging` (`SNP`, `A1`, `A2`,
#'   `Neighbours`, `Tagging` = \eqn{t_j}, `Weight`, `MAF`, `Exp_Heritability`,
#'   then one \eqn{v_{jc}} column per category), `ssums` (category x category:
#'   expected heritability that category `q2` contributes to SNPs of category
#'   `q`), `counts` (panel SNPs per category), `M` (panel SNPs) and `weights`.
#' @export
ldak.tagging <- function(eigen.files, snp.info, annot = NULL, alpha = -0.25, weights = c("thin", "none"),
                         thin.r2 = 0.98, n.ref = NULL, regression.snps = NULL, chunk = 2000) {
  snp.info <- as.data.frame(snp.info)
  M <- nrow(snp.info)
  A <- matrix(1, M, 1, dimnames = list(NULL, "Base"))
  if (!is.null(annot)) {
    annot <- as.data.frame(annot)
    Aa <- as.matrix(annot[match(snp.info$SNP, annot$SNP), -1, drop = FALSE])
    Aa[is.na(Aa)] <- 0
    A <- cbind(A, Aa)
  }
  Q <- ncol(A)
  reg <- if (is.null(regression.snps)) rep(TRUE, M) else snp.info$SNP %in% regression.snps
  f <- snp.info$Freq
  het <- 2 * f * (1 - f)
  if (is.numeric(weights)) { w <- weights; weights <- "user" } else { weights <- match.arg(weights); w <- rep(1, M) }

  blocks <- lapply(eigen.files, function(fl) { e <- read_eigen_bin(fl); e$U %*% diag(sqrt(pmax(e$lambda, 0)), e$k) })
  sizes <- vapply(blocks, nrow, 1L)
  if (sum(sizes) != M) stop("eigen blocks have ", sum(sizes), " SNPs but snp.info has ", M)
  ends <- cumsum(sizes); starts <- ends - sizes + 1
  bias <- if (is.null(n.ref)) 0 else 1 / (n.ref - 2)

  # LDAK-Thin weights
  if (weights == "thin") {
    for (b in seq_along(blocks)) {
      X <- blocks[[b]]; d <- rowSums(X^2); m <- nrow(X); keep <- rep(TRUE, m)
      for (s in seq(1, m, by = chunk)) {
        idx <- s:min(m, s + chunk - 1)
        r2 <- .block_r2(X, d, idx)
        for (ii in seq_along(idx)) {
          j <- idx[ii]
          if (j < m && keep[j] && any(keep[(j + 1):m] & r2[ii, (j + 1):m] >= thin.r2)) keep[j] <- FALSE
        }
      }
      w[starts[b]:ends[b]] <- as.numeric(keep)
    }
  }
  exps <- w * het^(1 + alpha)
  exps[!is.finite(exps) | het == 0] <- 0
  P <- A * exps                      # q_kc = a_kc w_k [2f(1-f)]^(1+alpha)

  V <- matrix(0, M, Q, dimnames = list(NULL, colnames(A)))
  tg <- numeric(M)
  for (b in seq_along(blocks)) {
    X <- blocks[[b]]; d <- rowSums(X^2); m <- nrow(X); rows <- starts[b]:ends[b]
    Pb <- P[rows, , drop = FALSE]; regb <- reg[rows]
    for (s in seq(1, m, by = chunk)) {
      idx <- s:min(m, s + chunk - 1)
      r2 <- .block_r2(X, d, idx)
      r2 <- r2 - (1 - r2) * bias
      r2[cbind(seq_along(idx), idx)] <- 1
      V[rows[idx], ] <- r2 %*% Pb
      tg[rows[idx]] <- as.vector(r2 %*% as.numeric(regb))
    }
  }

  member <- A != 0
  ssums <- crossprod(member * 1, P)  # [q, q2] = sum over SNPs in q of q_k,q2
  dimnames(ssums) <- list(colnames(A), colnames(A))
  tagging <- data.frame(SNP = snp.info$SNP, A1 = snp.info$A1, A2 = snp.info$A2,
                        Neighbours = rep(sizes, sizes), Tagging = tg, Weight = w, MAF = pmin(f, 1 - f),
                        Exp_Heritability = exps, V, check.names = FALSE)[reg, ]
  rownames(tagging) <- NULL
  structure(list(tagging = tagging, ssums = ssums, counts = colSums(member), M = M, M.reg = sum(reg),
                 alpha = alpha, weights = weights), class = "ldak_tagging")
}

#' Read an LDAK tagging file (output of `ldak --calc-tagging`)
#' @param file Path to the `.tagging` file.
#' @return An `ldak_tagging` list as from [ldak.tagging()]; the Base category
#'   (last in LDAK files) is moved first.
#' @export
read_ldak_tagging <- function(file) {
  lines <- readLines(file)
  hdr <- strsplit(lines[1], " +")[[1]]
  cats <- hdr[-(1:9)]; Q <- length(cats)
  nrow <- length(lines) - 2 - Q
  tab <- utils::read.table(text = lines[2:(nrow + 1)], col.names = hdr, check.names = FALSE, stringsAsFactors = FALSE,
                           colClasses = c("character", "character", "character", rep("numeric", 6 + Q)))
  ss <- t(sapply(lines[(nrow + 2):(nrow + 1 + Q)], function(x) as.numeric(utils::tail(strsplit(x, " +")[[1]], Q)), USE.NAMES = FALSE))
  ssums <- t(ss)  # file row q2 lists ssums[q][q2] over q
  last <- strsplit(lines[nrow + 2 + Q], " +")[[1]]
  counts <- as.numeric(utils::tail(last, Q)); Mref <- as.numeric(last[3])
  names(tab)[names(tab) == "Predictor"] <- "SNP"
  ord <- if (Q > 1) c(Q, seq_len(Q - 1)) else 1
  if (Q > 1) cats[Q] <- "Base"
  V <- as.matrix(tab[, 9 + seq_len(Q), drop = FALSE])[, ord, drop = FALSE]
  colnames(V) <- cats[ord]
  ssums <- ssums[ord, ord, drop = FALSE]; dimnames(ssums) <- list(cats[ord], cats[ord])
  tagging <- data.frame(tab[, c("SNP", "A1", "A2", "Neighbours", "Tagging", "Weight", "MAF", "Exp_Heritability")], V, check.names = FALSE)
  structure(list(tagging = tagging, ssums = ssums, counts = stats::setNames(counts[ord], cats[ord]), M = Mref, M.reg = nrow(tab),
                 alpha = NA, weights = NA), class = "ldak_tagging")
}

#' Write an `ldak_tagging` object in LDAK's tagging-file format
#'
#' The file can be used with `ldak --sum-hers` / `--sum-cors`. The Base
#' category is written last, as LDAK expects.
#' @param tag An `ldak_tagging` object.
#' @param file Output path.
#' @export
write_ldak_tagging <- function(tag, file) {
  Q <- ncol(tag$ssums); cats <- colnames(tag$ssums)
  ord <- if (Q > 1) c(seq_len(Q)[-1], 1) else 1
  tg <- tag$tagging
  V <- as.matrix(tg[, 8 + seq_len(Q), drop = FALSE])[, ord, drop = FALSE]
  lab <- if (Q > 1) c(paste0("Annotation", seq_len(Q - 1)), "Base") else "Base"
  out <- c(paste("Predictor A1 A2 Neighbours Tagging Weight MAF Categories Exp_Heritability", paste(lab, collapse = " ")),
           paste(tg$SNP, tg$A1, tg$A2, as.integer(tg$Neighbours), sprintf("%.3f", tg$Tagging), sprintf("%.4f", tg$Weight),
                 sprintf("%.6f", tg$MAF), -1, sprintf("%.4f", tg$Exp_Heritability),
                 apply(V, 1, function(v) paste(sprintf("%.4f", v), collapse = " "))))
  ss <- tag$ssums[ord, ord, drop = FALSE]
  for (q2 in seq_len(Q)) {
    what <- if (ord[q2] == 1) "the Base" else paste("Annotation", q2)
    out <- c(out, paste("The relative contribution of", what, "to each category", paste(sprintf("%.4f", ss[, q2]), collapse = " ")))
  }
  out <- c(out, paste("There are", tag$M, "reference", tag$M.reg, "regression", tag$M, "heritability predictors",
                      paste(as.integer(tag$counts[ord]), collapse = " ")))
  writeLines(out, file)
  invisible(file)
}

#' SumHer: heritability and enrichment from summary statistics (LDAK model)
#'
#' Maximises SumHer's approximate likelihood: each \eqn{S_j = Z_j^2} is a
#' scaled \eqn{\chi^2_1} with mean
#' \deqn{E[S_j] = 1 + n_j a + n_j \sum_c \theta_c v_{jc} ,}
#' and SNP \eqn{j} is down-weighted by its tagging \eqn{t_j} (over-counting),
#' \deqn{\ell = -\tfrac12 \sum_j \{ S_j / E_j + \log(S_j E_j) \} / t_j .}
#' Newton-Raphson with step halving, as in LDAK `--sum-hers`
#' (`--chisq-solver YES`). The heritability of category \eqn{c} is
#' \eqn{h^2_c = \theta_c \cdot \mathrm{ssums}_{cc}} and the heritability carried
#' by the SNPs of category \eqn{c} (including overlaps) is
#' \eqn{\sum_{c'} \theta_{c'} \mathrm{ssums}_{cc'}}.
#'
#' @param gwas `data.frame` with `SNP`, `Zscore`, `N` (aligned to the panel).
#' @param tag An `ldak_tagging` object ([ldak.tagging()] or [read_ldak_tagging()]).
#' @param intercept Estimate the additive confounding term \eqn{a}
#'   (LDAK `--intercept`; LDAK's default for `--sum-hers` is `FALSE`).
#' @param z.max Truncation of \eqn{|Z|}; `NULL` = \eqn{\sqrt{\max(80, 0.001 N_{max})}}, `Inf` = none (LDAK default).
#' @param se `"information"` (inverse observed information, LDAK default) or
#'   `"jackknife"` (delete-one-block over `nblock` contiguous blocks).
#' @param nblock Blocks for the jackknife.
#' @param tol,maxiter Convergence tolerance on the log-likelihood and iteration limit.
#' @return A list with `h2` (total, SE), `Heritability` (per category: `h2_c`),
#'   `Category` (heritability carried by the SNPs of each category),
#'   `Enrichment` (share, expected share = proportion of SNPs, enrichment),
#'   `intercept` (\eqn{1 + \bar n a}), `loglik`, `M.reg`.
#' @export
sumher <- function(gwas, tag, intercept = FALSE, z.max = NULL, se = c("information", "jackknife"), nblock = 200,
                   tol = 1e-3, maxiter = 100) {
  se <- match.arg(se)
  gwas <- as.data.frame(gwas); tg <- tag$tagging
  common <- tg$SNP[tg$SNP %in% gwas$SNP]
  if (length(common) < 0.8 * nrow(tg)) warning("GWAS covers ", round(100 * length(common) / nrow(tg)), "% of the regression SNPs; consider rebuilding the tagging with regression.snps = GWAS SNPs")
  tg <- tg[match(common, tg$SNP), ]; g <- gwas[match(common, gwas$SNP), ]
  Q <- ncol(tag$ssums); ss <- tag$ssums
  V <- as.matrix(tg[, 8 + seq_len(Q), drop = FALSE])
  if (is.null(z.max)) z.max <- .zmax_default(g$N)
  S <- pmin(g$Zscore^2, z.max^2); S[S == 0] <- 1e-6
  n <- g$N
  stag <- pmax(tg$Tagging, 1)
  scale <- sum(n / stag) / sum(1 / stag)
  Tm <- cbind(n / scale * V, if (intercept) n / scale)
  K <- ncol(Tm)
  .solve <- function(rows, theta) {
    S <- S[rows]; Tm <- Tm[rows, , drop = FALSE]; iw <- 1 / stag[rows]
    .lik <- function(E) -0.5 * sum((S / E + log(S * E)) * iw)
    .E <- function(th) pmax(1 + as.vector(Tm %*% th), 1e-6)
    E <- .E(theta); like <- .lik(E)
    for (it in seq_len(maxiter)) {
      g1 <- crossprod(Tm, 0.5 * (S - E) * iw / E^2)
      AI <- crossprod(Tm, Tm * ((S - 0.5 * E) * iw / E^3))
      step <- as.vector(solve(AI, g1))
      relax <- 1; moved <- FALSE
      while (relax > 1e-3) {
        th2 <- theta + relax * step; E2 <- .E(th2); l2 <- .lik(E2)
        if (l2 > like - tol) { moved <- TRUE; break }
        relax <- relax / 2
      }
      if (!moved) break
      diff <- l2 - like; theta <- th2; E <- E2; like <- l2
      if (abs(diff) < tol) break
    }
    AI <- crossprod(Tm, Tm * ((S - 0.5 * E) * iw / E^3))
    list(theta = theta, like = like, AI = AI)
  }
  fit <- .solve(seq_along(S), rep(0, K))
  th <- fit$theta
  # transforms: per-category h2, category-carried h2, total
  Jh <- cbind(diag(diag(ss) / scale, Q), matrix(0, Q, K - Q))
  Jc <- cbind(ss / scale, matrix(0, Q, K - Q))
  Jt <- colSums(Jh)
  est <- function(t) c(as.vector(Jh %*% t), as.vector(Jc %*% t), sum(Jt * t), if (intercept) t[K])
  e0 <- est(th)
  if (se == "information") {
    Vt <- solve(fit$AI)
    Jall <- rbind(Jh, Jc, Jt, if (intercept) c(rep(0, K - 1), 1))
    Vall <- Jall %*% Vt %*% t(Jall)
    sd0 <- sqrt(pmax(diag(Vall), 0))
    tot <- e0[2 * Q + 1]
    share.sd <- sapply(seq_len(Q), function(q) {
      gr <- (Jc[q, ] - e0[Q + q] / tot * Jt) / tot
      sqrt(max(as.numeric(t(gr) %*% Vt %*% gr), 0))
    })
  } else {
    IndexList <- .block_index(length(S), nblock); B <- length(IndexList)
    reps <- t(sapply(IndexList, function(idx) est(.solve(setdiff(seq_along(S), idx), th)$theta)))
    .jk <- function(v) sqrt((B - 1) / B * sum((v - mean(v))^2))
    sd0 <- apply(reps, 2, .jk)
    tot.r <- reps[, 2 * Q + 1]
    share.sd <- sapply(seq_len(Q), function(q) .jk(reps[, Q + q] / tot.r))
  }
  cats <- colnames(ss); tot <- e0[2 * Q + 1]
  expct <- tag$counts / tag$M
  list(h2 = data.frame(h2 = tot, se = sd0[2 * Q + 1]),
       Heritability = data.frame(Category = cats, h2 = e0[1:Q], se = sd0[1:Q]),
       Category = data.frame(Category = cats, h2 = e0[Q + 1:Q], se = sd0[Q + 1:Q]),
       Enrichment = data.frame(Category = cats, Share = e0[Q + 1:Q] / tot, Share.se = share.sd, Expected = expct,
                               Enrichment = e0[Q + 1:Q] / tot / expct, Enrichment.se = share.sd / expct),
       intercept = if (intercept) c(estimate = 1 + e0[length(e0)], se = sd0[length(sd0)]) else c(estimate = 1, se = 0),
       loglik = fit$like, z.max = z.max, M.reg = length(S), se.method = se)
}

#' Convert an `ldak_tagging` object to an LD-score table for the LDSC functions
#'
#' Under the LDAK model \eqn{E[\chi^2_j] = 1 + N h^2 v_j / \sum_k q_k}, so the
#' base tagging \eqn{v_j} plays the role of the LD score and
#' \eqn{\sum_k q_k = \mathrm{ssums}_{Base,Base}} the role of \eqn{M}. Use the
#' returned table with `M = attr(x, "M")` in [ldsc.univ()], [ldsc.bicov()] or
#' [ldsc.rg()] (e.g. genetic correlation under the LDAK-Thin model);
#' `wLD` = \eqn{t_j} gives the over-counting weight.
#' @param tag An `ldak_tagging` object.
#' @export
ldak.as.ldsc <- function(tag) {
  tg <- tag$tagging
  out <- data.frame(SNP = tg$SNP, LDSC = tg$Base, wLD = tg$Tagging)
  attr(out, "M") <- tag$ssums["Base", "Base"]
  out
}
