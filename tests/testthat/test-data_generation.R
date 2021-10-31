test_that("Data generation under default settings.", {
  withr::local_seed(101)
  data <- MCC::GenData()
  expect(is(data, "data.frame"), failure_message = "Data.frame not produced.")
})


test_that("Data generation with covariates but without betas.", {
  withr::local_seed(101)
  data0 <- MCC::GenData()
  
  n <- 100
  covariates <- stats::rnorm(n)
  withr::local_seed(101)
  data1 <- MCC::GenData(covariates = covariates)
  expect_identical(data0$time, data1$time)
  expect_identical(data0$status, data1$status)
})


test_that("Data generation with covariates and betas.", {
  withr::local_seed(101)
  n <- 100
  base_death_rate <- 1.0
  base_event_rate <- 1.0
  beta_death <- 1.0
  beta_event <- 1.0
  covariates <- stats::rnorm(n)
  data <- MCC::GenData(
    base_death_rate = base_death_rate,
    base_event_rate = base_event_rate,
    beta_death = beta_death,
    beta_event = beta_event,
    covariates = covariates
  )
  expect_gt(length(unique(data$true_death_rate)), 1)
  expect_gt(length(unique(data$true_event_rate)), 1)
})