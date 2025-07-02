#' @export
ldsc_bootstrap_wrapper <- function(ZMatrix1, NMatrix1, LDSC1,
                                   GCovInit, ECovInit,
                                   sampling.time = 200,
                                   zsquare_thresh = 50000,
                                   cov_thresh = 10000, n_threads = NULL) {

  # Convert to matrix
  Zmat <- as.matrix(as.data.frame(ZMatrix1[,-1]))
  Nmat <- as.matrix(as.data.frame(NMatrix1[,-1]))
  LDSC_vec <- LDSC1$LDSC
  Block_vec <- LDSC1$Block
  M <- nrow(Zmat)

  n_threads_cpp <- if (is.null(n_threads)) -1 else as.integer(n_threads)
  # Call C++ function
  res <- ldsc_bootstrap_cpp(
               Zmat, Nmat, LDSC_vec, Block_vec,
               as.matrix(GCovInit), as.matrix(ECovInit),
               as.integer(sampling.time),
               as.numeric(zsquare_thresh),
               as.numeric(cov_thresh),
               as.integer(M), n_threads_cpp)

  # Convert to array
  p <- ncol(Zmat)
  GCovEstt <- res$GCovEstt
  ECovEstt <- res$ECovEstt
  nsim <- length(GCovEstt)

  # Convert list of matrices → array
  GCovArray <- array(unlist(GCovEstt), dim = c(p, p, nsim))
  ECovArray <- array(unlist(ECovEstt), dim = c(p, p, nsim))
  GCorArray <- apply(GCovArray, 3, cov2cor_ldscr)
  GCorArray <- array(GCorArray, dim = c(p, p, nsim))

  # Compute standard errors
  GCovSE <- apply(GCovArray, c(1, 2), sd)
  ECovSE <- apply(ECovArray, c(1, 2), sd)
  GCorSE <- apply(GCorArray, c(1, 2), sd)

  return(list(
    GCovSE = GCovSE,
    ECovSE = ECovSE,
    GCorSE = GCorSE
  ))
}
