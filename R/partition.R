# Heritability partition shared by ldsc.strat() and sumher().
# tau: coefficients (Q); Cmat: Q x Q with h2_c = (Cmat tau)_c the heritability of the SNPs of category c;
# tvec: total h2 = sum(tvec * tau); Mc: SNPs per category; M: all SNPs; binary: which columns are 0/1.
# SEs from replicate coefficients (reps, R x Q, with se_fun) or from a covariance matrix V of tau (delta method).
.partition_table <- function(cats, tau, Cmat, tvec, Mc, M, binary, reps = NULL, se_fun = NULL, V = NULL) {
  pieces <- function(t) {
    h2c <- as.vector(Cmat %*% t); tot <- sum(tvec * t)
    rest <- M - Mc
    diff <- ifelse(rest > 0, h2c / Mc - (tot - h2c) / rest, NA_real_)
    list(h2c = h2c, prop = h2c / tot, diff = diff, tot = tot)
  }
  p0 <- pieces(tau)
  if (!is.null(reps)) {
    rr <- lapply(seq_len(nrow(reps)), function(r) pieces(reps[r, ]))
    sd_of <- function(nm) apply(do.call(rbind, lapply(rr, `[[`, nm)), 2, se_fun)
    tau.se <- apply(reps, 2, se_fun); h2.se <- sd_of("h2c"); prop.se <- sd_of("prop"); diff.se <- sd_of("diff")
  } else if (!is.null(V)) {
    q <- function(g) sqrt(max(as.numeric(t(g) %*% V %*% g), 0))
    tau.se <- sqrt(pmax(diag(V), 0))
    h2.se <- apply(Cmat, 1, q)
    prop.se <- sapply(seq_along(cats), function(c) q((Cmat[c, ] - p0$prop[c] * tvec) / p0$tot))
    diff.se <- sapply(seq_along(cats), function(c) {
      if (M - Mc[c] <= 0) return(NA_real_)
      q(Cmat[c, ] / Mc[c] - (tvec - Cmat[c, ]) / (M - Mc[c]))
    })
  } else {
    tau.se <- h2.se <- prop.se <- diff.se <- rep(NA_real_, length(cats))
  }
  ps <- Mc / M
  tz <- tau / tau.se
  out <- data.frame(Category = cats, Prop_SNPs = ps, h2 = p0$h2c, h2.se = h2.se, Prop_h2 = p0$prop, Prop_h2.se = prop.se,
                    Enrichment = p0$prop / ps, Enrichment.se = prop.se / ps,
                    Enrichment_p = 2 * stats::pnorm(-abs(p0$diff / diff.se)),
                    tau = tau, tau.se = tau.se, tau_z = tz, tau_p = 2 * stats::pnorm(-abs(tz)))
  nb <- !binary
  out[nb, c("Prop_SNPs", "h2", "h2.se", "Prop_h2", "Prop_h2.se", "Enrichment", "Enrichment.se", "Enrichment_p")] <- NA
  rownames(out) <- NULL
  out
}
