suppressPackageStartupMessages({
  library(dplyr)
})

#' Expect a ggplot or cowplot gtable plot object.
#'
#' @param p Plot object.
#' @noRd
.ExpectUniPlot <- function(p) {
  is_plot <- inherits(p, "ggplot") || inherits(p, "gtable")
  testthat::expect_true(is_plot)
}


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

  q_with <- PlotMCFs(data = data, strata_name = "strata", show_nar = FALSE)
  q_without <- PlotMCFs(data = data, show_nar = FALSE)

  expect_error(show(q_with), NA)
  expect_error(show(q_without), NA)
  expect_s3_class(q_with, "ggplot")
  expect_s3_class(q_without, "ggplot")
})

test_that("PlotMCFs default returns combined MCF + NAR layout.", {
  withr::local_seed(107)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(covariates = covariates, tau = 4)
  p <- PlotMCFs(data = data, tau = 4)
  .ExpectUniPlot(p)
  expect_s3_class(p, "ggplot")
  expect_error(show(p), NA)
})

test_that("PlotMCFs show_nar = FALSE returns MCF-only ggplot.", {
  withr::local_seed(108)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(covariates = covariates, tau = 4)
  p <- PlotMCFs(data = data, tau = 4, show_nar = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "top")
})

test_that("PlotMCFs which_arm plots one arm with combined layout.", {
  withr::local_seed(109)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(covariates = covariates, tau = 4)
  p <- PlotMCFs(
    data = data,
    which_arm = 0,
    tau = 4,
    x_breaks = 0:4,
    color_lab = "Control"
  )
  .ExpectUniPlot(p)
  expect_error(show(p), NA)
})

test_that("PlotMCFs show_auc runs for one- and two-arm data.", {
  withr::local_seed(110)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(covariates = covariates, tau = 4)
  expect_error(show(PlotMCFs(
    data = data,
    tau = 4,
    show_nar = FALSE,
    show_auc = TRUE
  )), NA)
  data_one <- GenData(n = 50, tau = 4)
  expect_error(show(PlotMCFs(
    data = data_one,
    tau = 4,
    show_nar = FALSE,
    show_auc = TRUE,
    color_lab = "Single Arm"
  )), NA)
})

test_that("PlotMCFs show_auc with show_nar runs without error.", {
  withr::local_seed(111)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(covariates = covariates, tau = 4)
  p <- PlotMCFs(
    data = data,
    tau = 4,
    x_breaks = 0:4,
    show_nar = TRUE,
    show_auc = TRUE
  )
  .ExpectUniPlot(p)
  expect_error(show(p), NA)
})

# -----------------------------------------------------------------------------
# Smoke tests for deprecated plot wrappers (no crash)
# -----------------------------------------------------------------------------

test_that("PlotOneSampleMCF runs without error.", {
  withr::local_seed(102)
  data <- MCC::GenData(n = 30, tau = 2)
  data$arm <- 1
  p <- PlotOneSampleMCF(data = data, tau = 2)
  expect_error(show(p), NA)
  expect_s3_class(p, "ggplot")
})

test_that("PlotOneSampleAUMCF runs without error.", {
  withr::local_seed(103)
  data <- MCC::GenData(n = 30, tau = 2)
  data$arm <- 1
  p <- PlotOneSampleAUMCF(data = data, tau = 2)
  expect_error(show(p), NA)
  expect_s3_class(p, "ggplot")
})

test_that("PlotOneSampleNAR runs without error.", {
  withr::local_seed(104)
  data <- MCC::GenData(n = 30, tau = 2)
  data$arm <- 1
  x_breaks <- seq(0, 2, by = 0.5)
  p <- PlotOneSampleNAR(data = data, x_breaks = x_breaks)
  expect_error(show(p), NA)
  expect_s3_class(p, "ggplot")
})

test_that("PlotNARs runs without error.", {
  withr::local_seed(105)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- MCC::GenData(covariates = covariates, tau = 2)
  x_breaks <- seq(0, 2, by = 0.5)
  p <- PlotNARs(data = data, x_breaks = x_breaks)
  expect_error(show(p), NA)
  expect_s3_class(p, "ggplot")
})

test_that("PlotAUMCF runs without error.", {
  withr::local_seed(106)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- MCC::GenData(covariates = covariates, tau = 2)
  p <- PlotAUMCF(data = data, which_arm = 1, tau = 2)
  expect_error(show(p), NA)
  expect_s3_class(p, "ggplot")
})
