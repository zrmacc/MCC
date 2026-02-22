suppressPackageStartupMessages({
  library(dplyr)
})

test_that("Plotting with and without stratification.", {
  withr::local_seed(101)
  covariates <- data.frame(
    arm = rep(c(0, 1), each = 50),
    strata = rbinom(n = 100, size = 2, prob = 0.25)
  )
  
  data <- MCC::GenData(
    covariates = covariates,
    beta_event = c(log(0.8), log(1.2))
  )
  
  # Plotting with and without strata.
  q_with <- PlotMCFs(data = data, strata_name = "strata")
  q_without <- PlotMCFs(data = data)
  
  expect_error(show(q_with), NA)
  expect_error(show(q_without), NA)
})

# -----------------------------------------------------------------------------
# Smoke tests for other plot functions (no crash)
# -----------------------------------------------------------------------------

test_that("PlotOneSampleMCF runs without error.", {
  withr::local_seed(102)
  data <- MCC::GenData(n = 30, tau = 2)
  data$arm <- 1
  p <- PlotOneSampleMCF(data = data, tau = 2)
  expect_error(show(p), NA)
})

test_that("PlotOneSampleAUMCF runs without error.", {
  withr::local_seed(103)
  data <- MCC::GenData(n = 30, tau = 2)
  data$arm <- 1
  p <- PlotOneSampleAUMCF(data = data, tau = 2)
  expect_error(show(p), NA)
})

test_that("PlotOneSampleNAR runs without error.", {
  withr::local_seed(104)
  data <- MCC::GenData(n = 30, tau = 2)
  data$arm <- 1
  x_breaks <- seq(0, 2, by = 0.5)
  p <- PlotOneSampleNAR(data = data, x_breaks = x_breaks)
  expect_error(show(p), NA)
})

test_that("PlotNARs runs without error.", {
  withr::local_seed(105)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- MCC::GenData(covariates = covariates, tau = 2)
  x_breaks <- seq(0, 2, by = 0.5)
  p <- PlotNARs(data = data, x_breaks = x_breaks)
  expect_error(show(p), NA)
})

test_that("PlotAUMCF runs without error.", {
  withr::local_seed(106)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- MCC::GenData(covariates = covariates, tau = 2)
  p <- PlotAUMCF(data = data, which_arm = 1, tau = 2)
  expect_error(show(p), NA)
})

