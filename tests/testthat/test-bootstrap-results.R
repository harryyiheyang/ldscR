bootstrap_test_data <- function() {
  M <- 120L
  index <- seq_len(M)
  ldsc <- data.frame(
    SNP = paste0("rs", index),
    LDSC = seq(1, 3, length.out = M)
  )
  N <- 10000
  l <- ldsc$LDSC * N / M

  list(
    ldsc = ldsc,
    gwas1 = data.frame(
      SNP = ldsc$SNP,
      Zscore = sqrt(1.05 + 0.01 * l + 0.05 * sin(index / 4)),
      N = N,
      A1 = "A",
      A2 = "G"
    ),
    gwas2 = data.frame(
      SNP = ldsc$SNP,
      Zscore = sqrt(1.10 + 0.012 * l + 0.04 * cos(index / 5)),
      N = N,
      A1 = "A",
      A2 = "G"
    )
  )
}

test_that("ldsc.univ exposes the existing bootstrap results", {
  inputs <- bootstrap_test_data()
  sampling.time <- 24L

  set.seed(20260823)
  invisible(capture.output(
    fit <- ldsc.univ(
      inputs$gwas1,
      inputs$ldsc,
      sampling.time = sampling.time,
      nblock = 8,
      two.step = NULL,
      z.max = Inf
    )
  ))

  expect_true(is.data.frame(fit))
  expect_identical(
    names(fit),
    c("intercept", "intercept.se", "h2", "h2.se", "M")
  )
  expect_length(fit$bootstrap$intercept, sampling.time)
  expect_length(fit$bootstrap$h2, sampling.time)
  expect_equal(sd(fit$bootstrap$intercept), fit$intercept.se)
  expect_equal(sd(fit$bootstrap$h2), fit$h2.se)
  expect_identical(dim(fit$bootstrap$block.counts), c(sampling.time, 8L))
  expect_true(all(rowSums(fit$bootstrap$block.counts) == 8L))

  # Point-estimate fixtures recorded from commit 1a61a3c7514e18e9be6f7ba51b6a3b5cc243d9c0.
  expect_equal(fit$intercept, 1.07480674790884, tolerance = 1e-12)
  expect_equal(fit$h2, 0.00985315689304359, tolerance = 1e-12)
  expect_equal(fit$intercept.se, 0.0320863349100051, tolerance = 1e-12)
  expect_equal(fit$h2.se, 0.000211270195393837, tolerance = 1e-12)
})

test_that("ldsc.bicov exposes the existing bootstrap results", {
  inputs <- bootstrap_test_data()
  sampling.time <- 24L

  set.seed(20260823)
  invisible(capture.output(
    fit <- ldsc.bicov(
      inputs$gwas1,
      inputs$gwas2,
      LDSC = inputs$ldsc,
      sampling.time = sampling.time,
      nblock = 8,
      two.step = NULL,
      z.max = Inf
    )
  ))

  expect_true(is.data.frame(fit))
  expect_identical(names(fit), c("ecov", "ecov.se", "gcov", "gcov.se", "h21", "h21.se", "h22", "h22.se", "rg", "rg.se", "M"))
  expect_length(fit$bootstrap$ecov, sampling.time)
  expect_length(fit$bootstrap$gcov, sampling.time)
  expect_equal(sd(fit$bootstrap$ecov), fit$ecov.se)
  expect_equal(sd(fit$bootstrap$gcov), fit$gcov.se)
  expect_identical(dim(fit$bootstrap$block.counts), c(sampling.time, 8L))
  expect_true(all(rowSums(fit$bootstrap$block.counts) == 8L))

  # Fixtures re-recorded: Var(z1 z2) = v1 v2 + mu^2 weights, h21 and h22 estimated inside ldsc.bicov.
  expect_equal(fit$ecov, 1.09194136966973, tolerance = 1e-12)
  expect_equal(fit$gcov, 0.0108557833818878, tolerance = 1e-12)
  expect_equal(fit$ecov.se, 0.0187300476150324, tolerance = 1e-12)
  expect_equal(fit$gcov.se, 0.000113372728088476, tolerance = 1e-12)
})
