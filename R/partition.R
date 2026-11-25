# Heritability partition shared by ldsc.strat() and sumher().
# tau: coefficients (Q); tvec: per-column scale with the additive h2 of column c = tau_c * tvec_c and total
# h2 = sum_c tau_c tvec_c (tvec = sum_j a_jc for S-LDSC, ssums_cc for SumHer); Cmat: Q x Q with
# (Cmat tau)_c = h2 of the SNPs of category c including overlapping columns (ldsc --overlap-annot);
# Mc: SNPs per category; M: all SNPs; binary: which columns are 0/1.
# SEs from replicate coefficients (reps, R x Q, with se_fun) or from a covariance matrix V of tau (delta method).
.partition_table <- function(cats, tau, Cmat, tvec, Mc, M, binary, reps = NULL, se_fun = NULL, V = NULL) {
  pieces <- function(t) {
    own <- t * tvec; tot <- sum(own)
    h2c <- as.vector(Cmat %*% t)
    rest <- M - Mc
    diff <- ifelse(rest > 0, h2c / Mc - (tot - h2c) / rest, NA_real_)
    list(own = own, prop_own = own / tot, h2c = h2c, prop = h2c / tot, diff = diff, tot = tot)
  }
  p0 <- pieces(tau)
  Q <- length(cats)
  if (!is.null(reps)) {
    rr <- lapply(seq_len(nrow(reps)), function(r) pieces(reps[r, ]))
    sd_of <- function(nm) apply(do.call(rbind, lapply(rr, `[[`, nm)), 2, se_fun)
    tau.se <- apply(reps, 2, se_fun)
    own.se <- sd_of("own"); pown.se <- sd_of("prop_own")
    h2.se <- sd_of("h2c"); prop.se <- sd_of("prop"); diff.se <- sd_of("diff")
  } else if (!is.null(V)) {
    q <- function(g) sqrt(max(as.numeric(t(g) %*% V %*% g), 0))
    tau.se <- sqrt(pmax(diag(V), 0))
    own.se <- tau.se * abs(tvec)
    pown.se <- sapply(seq_len(Q), function(c) { g <- -p0$prop_own[c] * tvec / p0$tot; g[c] <- g[c] + tvec[c] / p0$tot; q(g) })
    h2.se <- apply(Cmat, 1, q)
    prop.se <- sapply(seq_len(Q), function(c) q((Cmat[c, ] - p0$prop[c] * tvec) / p0$tot))
    diff.se <- sapply(seq_len(Q), function(c) {
      if (M - Mc[c] <= 0) return(NA_real_)
      q(Cmat[c, ] / Mc[c] - (tvec - Cmat[c, ]) / (M - Mc[c]))
    })
  } else {
    tau.se <- own.se <- pown.se <- h2.se <- prop.se <- diff.se <- rep(NA_real_, Q)
  }
  ps <- Mc / M
  tz <- tau / tau.se
  out <- data.frame(Category = cats, Prop_SNPs = ps,
                    h2 = p0$own, h2.se = own.se, Prop_h2 = p0$prop_own, Prop_h2.se = pown.se,
                    tau = tau, tau.se = tau.se, tau_z = tz, tau_p = 2 * stats::pnorm(-abs(tz)),
                    h2_overlap = p0$h2c, h2_overlap.se = h2.se, Prop_h2_overlap = p0$prop, Prop_h2_overlap.se = prop.se,
                    Enrichment_overlap = p0$prop / ps, Enrichment_overlap.se = prop.se / ps,
                    Enrichment_overlap_p = 2 * stats::pnorm(-abs(p0$diff / diff.se)))
  ov <- grep("overlap", names(out), value = TRUE)
  out[!binary, ov] <- NA
  rownames(out) <- NULL
  out
}
