#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @export
ldsc.strat = function(gwas, SLDSC, annot, sampling.time = 0, nblock = 500, intercetp = T, intercetp.lower = 0.9, intercetp.upper = 2) {
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
  t0 <- Sys.time() - t0
  print("Processing data"); print(t0)

  X = as.matrix(SLDSC[, -1])
  Annot = as.matrix(annot[, -1])
  annot_names = colnames(annot)[-1]
  Annot_cov = colSums(Annot)
  remove(SLDSC, annot)

  .wls_fit <- function(X, z, w, intercetp = T, int.low = 0.9, int.up = 2) {
    if (intercetp == T) {
      X <- cbind(1, X)
      XtX <- CppMatrix::matrixMultiply(t(X), X * w)
      Xty <- CppMatrix::matrixVectorMultiply(t(X), z * w)
      res <- as.vector(qr.solve(XtX, Xty))

      if (res[1] < int.low) {
        res[1] <- int.low
        res[-1] <- as.vector(qr.solve(XtX[-1, -1, drop = FALSE], Xty[-1] - int.low * XtX[-1, 1]))
      } else if (res[1] > int.up) {
        res[1] <- int.up
        res[-1] <- as.vector(qr.solve(XtX[-1, -1, drop = FALSE], Xty[-1] - int.up * XtX[-1, 1]))
      }
      res
    } else {
      XtX <- CppMatrix::matrixMultiply(t(X), X * w)
      Xty <- CppMatrix::matrixVectorMultiply(t(X), (z - 1) * w)
      c(1, as.vector(qr.solve(XtX, Xty)))
    }
  }

  t1 <- Sys.time()
  z <- gwas$Zscore^2
  L_mat <- X * gwas$N
  l_vec <- X[, 1] * gwas$N / M
  w <- 1 / pmax(1, X[, 1]) / (1 + 0.1 * l_vec)^2

  beta_ini <- .wls_fit(L_mat, z, w, intercetp = intercetp, int.low = intercetp.lower, int.up = intercetp.upper)
  tau_ini <- beta_ini[-1]
  h2.ini <- sum(tau_ini * Annot_cov)
  t1 <- Sys.time() - t1
  print("Initial Heritability Estimate"); print(t1)

  t2 <- Sys.time()
  h2_curr <- h2.ini
  for (i in 1:3) {
    w <- 1 / (1 + l_vec * h2_curr)^2 / pmax(1, X[, 1])
    beta_tmp <- .wls_fit(L_mat, z, w, intercetp = intercetp, int.low = intercetp.lower, int.up = intercetp.upper)
    intercept <- beta_tmp[1]
    tau <- beta_tmp[-1]
    h2_curr <- sum(tau * Annot_cov)
  }
  h2 <- h2_curr
  t2 <- Sys.time() - t2
  print("Heritability Estimation (Reweighted)"); print(t2)

  if (sampling.time > 0) {
    t3 <- Sys.time()

    K <- if(intercetp) ncol(L_mat) + 1 else ncol(L_mat)

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
    ends <- starts + block_lengths - 1

    XtXlist <- vector("list", B)
    Xtylist <- vector("list", B)

    for (b in seq_len(B)) {
      idx <- starts[b]:ends[b]
      Xi <- L_mat[idx, , drop = FALSE]
      wi <- w[idx]
      zi <- z[idx]

      if (intercetp) {
        Xi <- cbind(1, Xi)
        XtXlist[[b]] <- CppMatrix::matrixMultiply(t(Xi), Xi * wi)
        Xtylist[[b]] <- CppMatrix::matrixVectorMultiply(t(Xi), zi * wi)
      } else {
        XtXlist[[b]] <- CppMatrix::matrixMultiply(t(Xi), Xi * wi)
        Xtylist[[b]] <- CppMatrix::matrixVectorMultiply(t(Xi), (zi - 1) * wi)
      }
    }

    h2.vec <- numeric(sampling.time)
    intercept.vec <- numeric(sampling.time)
    tau.mat <- matrix(NA, nrow = sampling.time, ncol = length(tau))

    for (r in seq_len(sampling.time)) {
      draw <- sample.int(B, B, replace = TRUE)
      cnt  <- tabulate(draw, nbins = B)
      sel  <- which(cnt > 0L)

      XtXi_sum <- matrix(0, K, K)
      Xtyi_sum <- numeric(K)

      for (j in sel) {
        cj <- cnt[j]
        XtXi_sum <- XtXi_sum + cj * XtXlist[[j]]
        Xtyi_sum <- Xtyi_sum + cj * Xtylist[[j]]
      }

      res_beta <- as.vector(qr.solve(XtXi_sum, Xtyi_sum))

      if (intercetp) {
        if (res_beta[1] < intercetp.lower) {
          res_beta[1] <- intercetp.lower
          res_beta[-1] <- as.vector(qr.solve(XtXi_sum[-1, -1, drop = FALSE], Xtyi_sum[-1] - intercetp.lower * XtXi_sum[-1, 1]))
        } else if (res_beta[1] > intercetp.upper) {
          res_beta[1] <- intercetp.upper
          res_beta[-1] <- as.vector(qr.solve(XtXi_sum[-1, -1, drop = FALSE], Xtyi_sum[-1] - intercetp.upper * XtXi_sum[-1, 1]))
        }
        intercept.vec[r] <- res_beta[1]
        tau_boot <- res_beta[-1]
      } else {
        intercept.vec[r] <- 1
        tau_boot <- res_beta
      }

      tau.mat[r, ] <- tau_boot
      h2.vec[r] <- sum(tau_boot * Annot_cov)
    }

    t3 <- Sys.time() - t3
    print("Standard Error Estimation"); print(t3)
    intercept.se <- sd(intercept.vec)
    h2.se <- sd(h2.vec)
    tau.se <- apply(tau.mat, 2, sd)
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
    Parameter = c("Intercept","Residual", annot_names[-1]),
    Estimate = c(intercept, tau),
    SE = c(intercept.se, tau.se)
  )

  list(H2_summary = H2_summary, Parameter_summary = Parameter_summary)
}
