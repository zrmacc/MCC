test_that("Test AUC function.", {
  
  # Case 1.
  times <- c(1, 2)
  values <- c(1, 2)
  tau <- 3
  observed <- AUC(times, values, tau)
  expected <- (tau - 1) + (tau - 2)
  expect_equal(observed, expected)
  
  # Case 2.
  times <- c(2, 4)
  values <- c(1, 3)
  tau <- 5
  observed <- AUC(times, values, tau)
  expected <- (tau - 2) + 2 * (tau - 4)
  expect_equal(observed, expected)
  
  # Case 3.
  times <- c(2, 4)
  values <- c(1, 3)
  tau <- 1
  observed <- AUC(times, values, tau)
  expected <- 0
  expect_equal(observed, expected)
  
  # Case 4.
  times <- c(1, 3)
  values <- c(2, 5)
  tau <- 3
  observed <- AUC(times, values, tau)
  expected <- 2 * (3 - 1) 
  expect_equal(observed, expected)
  
  # Case 5.
  times <- c(1, 4)
  values <- c(2, 5)
  tau <- 2
  observed <- AUC(times, values, tau)
  expected <- 2 * (2 - 1)
  expect_equal(observed, expected)
  
  # Case 6.
  times <- c(1, 2)
  values <- c(3, 5)
  tau <- 0
  observed <- AUC(times, values, tau)
  expected <- 0
  expect_equal(observed, expected)
  
  # Case 7.
  times <- c(0, 2)
  values <- c(1, 3)
  tau <- 3
  observed <- AUC(times, values, tau)
  expected <- 1 * (2 - 0) + 3 * (3 - 2)
  expect_equal(observed, expected)
  
  # Case 8.
  times <- c(1, 2, 3)
  values <- c(5, 5, 5)
  tau <- 4
  observed <- AUC(times, values, tau)
  expected <- 5 * (4 - 1)
  expect_equal(observed, expected)
})

# -----------------------------------------------------------------------------

test_that("AUC errors on non-monotone increasing values.", {
  times <- c(1, 2, 3)
  values <- c(1, 3, 2)
  tau <- 4
  expect_error(AUC(times, values, tau), "monotone increasing")
})

# -----------------------------------------------------------------------------
# VarAUC
# -----------------------------------------------------------------------------

test_that("VarAUC returns non-negative scalar variance.", {
  withr::local_seed(202)
  covar <- data.frame(arm = c(rep(1, 30), rep(0, 30)))
  data <- GenData(beta_event = log(0.5), covariates = covar, tau = 3)
  v <- VarAUC(data, tau = 2)
  expect_type(v, "double")
  expect_length(v, 1)
  expect_gte(v, 0)
  expect_true(is.finite(v))
})

test_that("VarAUC with return_psi = TRUE returns influence contributions.", {
  withr::local_seed(203)
  covar <- data.frame(arm = c(rep(1, 25), rep(0, 25)))
  data <- GenData(beta_event = log(0.5), covariates = covar, tau = 3)
  out <- VarAUC(data, tau = 2, return_psi = TRUE)
  expect_s3_class(out, "data.frame")
  expect_true("idx" %in% names(out))
  expect_true("psi" %in% names(out))
  expect_equal(nrow(out), length(unique(data$idx)))
})

test_that("VarAUC variance equals mean(psi^2) when return_psi = TRUE.", {
  withr::local_seed(204)
  covar <- data.frame(arm = c(rep(1, 20), rep(0, 20)))
  data <- GenData(beta_event = log(0.5), covariates = covar, tau = 3)
  v_direct <- VarAUC(data, tau = 2)
  out <- VarAUC(data, tau = 2, return_psi = TRUE)
  expect_equal(v_direct, mean(out$psi^2))
})
