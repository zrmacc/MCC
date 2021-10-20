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

