testthat::test_that("CalcAugComp returns within-arm centered sigma/gamma and xbar with correct shapes", {
  set.seed(1)
  
  n <- 40
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  psi <- rnorm(n)
  
  out <- CalcAugComp(X, psi)
  
  # Shapes.
  testthat::expect_equal(dim(out$sigma), c(p, p))
  testthat::expect_equal(length(out$gamma), p)
  testthat::expect_equal(dim(out$xbar), c(p, 1))
  
  # xbar equals column means.
  testthat::expect_equal(as.numeric(out$xbar), colMeans(X), tolerance = 1e-12)
  
  # Sigma calculation.
  resid_manual <- sweep(X, 2, colMeans(X), "-")
  sigma_manual <- crossprod(resid_manual) / (n^2)
  testthat::expect_equal(out$sigma, sigma_manual, tolerance = 1e-12)
  
  # Gamma calculation.
  gamma_manual <- crossprod(resid_manual, psi) / (n^2)
  testthat::expect_equal(as.numeric(out$gamma), as.numeric(gamma_manual), tolerance = 1e-12)
})


testthat::test_that("Augmentation correction uses imbalance in xbar (not zeroed by within-arm centering)", {
  set.seed(2)
  
  n0 <- 50
  n1 <- 60
  p <- 4
  
  # Force a known mean difference
  shift <- c(0.5, -1.0, 0.0, 2.0)
  
  X0 <- matrix(rnorm(n0 * p), n0, p)
  X1 <- matrix(rnorm(n1 * p), n1, p) + matrix(rep(shift, each = n1), n1, p)
  
  psi0 <- rnorm(n0)
  psi1 <- rnorm(n1)
  
  a0 <- CalcAugComp(X0, psi0)
  a1 <- CalcAugComp(X1, psi1)
  
  sigma <- a0$sigma + a1$sigma
  gamma <- a0$gamma + a1$gamma
  imbalance <- a1$xbar - a0$xbar
  
  # imbalance should be close to shift (sampling noise allowed)
  testthat::expect_equal(as.numeric(imbalance), shift, tolerance = 0.35)
  
  omega <- MASS::ginv(sigma) %*% gamma
  correction <- as.numeric(crossprod(omega, imbalance))
  
  # If imbalance is nonzero and omega nonzero, correction should not be identically 0
  testthat::expect_true(is.finite(correction))
  testthat::expect_true(abs(correction) > 1e-10)
})


testthat::test_that("Augmentation correction is invariant to adding a constant to all covariates", {
  set.seed(3)
  
  n0 <- 40
  n1 <- 45
  p <- 3
  
  X0 <- matrix(rnorm(n0 * p), n0, p)
  X1 <- matrix(rnorm(n1 * p), n1, p) + 0.4
  
  psi0 <- rnorm(n0)
  psi1 <- rnorm(n1)
  
  const <- rep(10, p)
  
  a0 <- CalcAugComp(X0, psi0)
  a1 <- CalcAugComp(X1, psi1)
  sigma <- a0$sigma + a1$sigma
  gamma <- a0$gamma + a1$gamma
  imbalance <- a1$xbar - a0$xbar
  omega <- MASS::ginv(sigma) %*% gamma
  corr1 <- as.numeric(crossprod(omega, imbalance))
  
  a0c <- CalcAugComp(sweep(X0, 2, const, "+"), psi0)
  a1c <- CalcAugComp(sweep(X1, 2, const, "+"), psi1)
  sigmac <- a0c$sigma + a1c$sigma
  gammac <- a0c$gamma + a1c$gamma
  imbalance_c <- a1c$xbar - a0c$xbar
  omegac <- MASS::ginv(sigmac) %*% gammac
  corr2 <- as.numeric(crossprod(omegac, imbalance_c))
  
  # Within-arm centering makes sigma/gamma unchanged by constant shift;
  # imbalance also unchanged; so correction invariant.
  testthat::expect_equal(sigmac, sigma, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(gammac), as.numeric(gamma), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(imbalance_c), as.numeric(imbalance), tolerance = 1e-12)
  testthat::expect_equal(corr2, corr1, tolerance = 1e-10)
})


testthat::test_that("Augmentation correction matches a 'direct' computation using raw means and within-arm centered matrices", {
  set.seed(4)
  
  n0 <- 80
  n1 <- 75
  p <- 6
  
  X0 <- matrix(rnorm(n0 * p), n0, p)
  X1 <- matrix(rnorm(n1 * p), n1, p) + 0.2
  
  psi0 <- rnorm(n0)
  psi1 <- rnorm(n1)
  
  # Via CalcAugComp
  a0 <- CalcAugComp(X0, psi0)
  a1 <- CalcAugComp(X1, psi1)
  sigma <- a0$sigma + a1$sigma
  gamma <- a0$gamma + a1$gamma
  imbalance <- a1$xbar - a0$xbar
  omega <- MASS::ginv(sigma) %*% gamma
  corr_calc <- as.numeric(crossprod(omega, imbalance))
  
  # Manual "direct" computation
  r0 <- sweep(X0, 2, colMeans(X0), "-")
  r1 <- sweep(X1, 2, colMeans(X1), "-")
  
  sigma_direct <- crossprod(r0) / (n0^2) + crossprod(r1) / (n1^2)
  gamma_direct <- crossprod(r0, psi0) / (n0^2) + crossprod(r1, psi1) / (n1^2)
  imbalance_direct <- colMeans(X1) - colMeans(X0)
  omega_direct <- MASS::ginv(sigma_direct) %*% gamma_direct
  corr_direct <- as.numeric(crossprod(omega_direct, imbalance_direct))
  
  testthat::expect_equal(sigma, sigma_direct, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(gamma), as.numeric(gamma_direct), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(imbalance), as.numeric(imbalance_direct), tolerance = 1e-12)
  testthat::expect_equal(corr_calc, corr_direct, tolerance = 1e-10)
})

