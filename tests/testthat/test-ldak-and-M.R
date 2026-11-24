test_that("M defaults to the panel size, so dropping GWAS SNPs does not rescale h2", {
  set.seed(5)
  M <- 4000
  ld <- data.frame(SNP = paste0("rs", 1:M), LDSC = rgamma(M, 4, 0.5))
  N <- 1e5; h2 <- 0.3
  z <- rnorm(M) * sqrt(1 + N * h2 * ld$LDSC / M)
  g <- data.frame(SNP = ld$SNP, Zscore = z, N = N)
  keep <- sort(sample.int(M, 0.7 * M))
  invisible(capture.output({
    f1 <- ldsc.univ(g, ld, two.step = NULL, z.max = Inf)
    f2 <- ldsc.univ(g[keep, ], ld, two.step = NULL, z.max = Inf)
  }))
  expect_equal(f2$M, M)
  expect_equal(f2$M.reg, length(keep))
  expect_lt(abs(f2$h2 / f1$h2 - 1), 0.25)
})

test_that("eigen blocks round-trip and LDAK taggings match a direct computation", {
  set.seed(7)
  m <- 60; n <- 500
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  for (j in 2:m) G[, j] <- ifelse(runif(n) < 0.7, G[, j - 1], G[, j])
  R <- cor(G); e <- eigen(R, symmetric = TRUE)
  f <- tempfile(fileext = ".eigen.bin")
  write_eigen_bin(f, e$vectors, e$values)
  eb <- read_eigen_bin(f)
  expect_equal(eb$U %*% diag(eb$lambda) %*% t(eb$U), R, tolerance = 1e-5)

  info <- data.frame(SNP = paste0("s", 1:m), A1 = "A", A2 = "C", Freq = colMeans(G) / 2)
  tag <- ldak.tagging(f, info, alpha = -0.25, weights = "none", n.ref = n)
  q <- (2 * info$Freq * (1 - info$Freq))^0.75
  r2 <- R^2; r2 <- r2 - (1 - r2) / (n - 2); diag(r2) <- 1
  expect_equal(tag$tagging$Base, as.vector(r2 %*% q), tolerance = 1e-4)
  expect_equal(tag$tagging$Tagging, rowSums(r2), tolerance = 1e-4)
  expect_equal(unname(tag$ssums[1, 1]), sum(q), tolerance = 1e-8)

  tf <- tempfile(fileext = ".tagging")
  write_ldak_tagging(tag, tf)
  back <- read_ldak_tagging(tf)
  expect_equal(back$tagging$Base, tag$tagging$Base, tolerance = 1e-3)
  expect_equal(unname(back$ssums), unname(tag$ssums), tolerance = 1e-3)
})

test_that("sumher recovers h2 under its own model", {
  set.seed(9)
  M <- 20000; N <- 2e5; h2 <- 0.25
  v <- rgamma(M, 3, 1); ssum <- M
  tag <- structure(list(
    tagging = data.frame(SNP = paste0("rs", 1:M), A1 = "A", A2 = "C", Neighbours = 1, Tagging = v, Weight = 1,
                         MAF = 0.3, Exp_Heritability = 1, Base = v),
    ssums = matrix(ssum, 1, 1, dimnames = list("Base", "Base")), counts = c(Base = M), M = M, M.reg = M), class = "ldak_tagging")
  E <- 1 + N * h2 * v / ssum
  g <- data.frame(SNP = tag$tagging$SNP, Zscore = rnorm(M) * sqrt(E), N = N)
  r <- sumher(g, tag, z.max = Inf)
  expect_lt(abs(r$h2$h2 - h2), 4 * r$h2$se)
})

test_that("ldsc.bicov estimates h2 internally and gives rg with a bootstrap SE", {
  set.seed(3)
  M <- 6000; N <- 5e4
  ld <- data.frame(SNP = paste0("rs", 1:M), LDSC = rgamma(M, 4, 0.5))
  s <- sqrt(1 + N * 0.3 * ld$LDSC / M)
  z1 <- rnorm(M) * s; z2 <- 0.6 * z1 + sqrt(1 - 0.36) * rnorm(M) * s
  g1 <- data.frame(SNP = ld$SNP, Zscore = z1, N = N); g2 <- data.frame(SNP = ld$SNP, Zscore = z2, N = N)
  set.seed(1)
  f <- ldsc.bicov(g1, g2, LDSC = ld, sampling.time = 50, nblock = 50)
  set.seed(1)
  u <- ldsc.univ(g1, ld, sampling.time = 50, nblock = 50)
  expect_equal(f$h21, u$h2, tolerance = 1e-10)
  expect_equal(f$h21.se, u$h2.se, tolerance = 1e-10)
  expect_true(is.finite(f$rg.se) && f$rg.se > 0)
  expect_equal(f$rg, f$gcov / sqrt(f$h21 * f$h22))
})
