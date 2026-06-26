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
  expect_gt(length(unique(data$death_rate)), 1)
  expect_gt(length(unique(data$event_rate)), 1)
})


test_that("K = 1 frailty_variance aliases death_frailty_variance.", {
  withr::local_seed(99)
  data1 <- MCC::GenData(frailty_variance = 1)
  withr::local_seed(99)
  data2 <- MCC::GenData(death_frailty_variance = 1)
  expect_equal(data1, data2)
})


test_that("Multivariate GenData produces valid long format.", {
  withr::local_seed(42)
  n <- 40
  covariates <- data.frame(arm = c(rep(1, n / 2), rep(0, n / 2)))
  data <- MCC::GenData(
    covariates = covariates,
    n_event_types = 2L,
    n = n
  )

  expect_true("event_type" %in% names(data))
  recurrent <- data$status == 1L
  terminal <- data$status %in% c(0L, 2L)
  expect_false(any(is.na(data$event_type[recurrent])))
  expect_true(all(is.na(data$event_type[terminal])))
  expect_equal(sort(unique(data$event_type[recurrent])), c(1, 2))
  terminal_counts <- table(data$idx[terminal])
  expect_true(all(as.integer(terminal_counts) == 1L))
})


test_that("Legacy frailty_variance warns when n_event_types > 1.", {
  expect_warning(
    MCC::GenData(n = 10, n_event_types = 2L, frailty_variance = 1),
    "death_frailty_variance and event_frailty_variance"
  )
})


test_that("Death frailty scales rates when event frailty is zero.", {
  withr::local_seed(11)
  data <- MCC::GenData(
    n = 50,
    n_event_types = 2L,
    death_frailty_variance = 2,
    event_frailty_variance = 0
  )
  subj <- unique(data[, c("idx", "death_frailty", "event_frailty")])
  expect_gt(length(unique(subj$death_frailty)), 1L)
  expect_true(all(subj$event_frailty == 1))
})


test_that("Event frailty varies when death frailty is zero.", {
  withr::local_seed(12)
  data <- MCC::GenData(
    n = 50,
    n_event_types = 2L,
    death_frailty_variance = 0,
    event_frailty_variance = 2
  )
  subj <- unique(data[, c("idx", "death_frailty", "event_frailty")])
  expect_true(all(subj$death_frailty == 1))
  expect_gt(length(unique(subj$event_frailty)), 1L)
})


test_that("CompareMvAUCs runs on multivariate GenData output.", {
  withr::local_seed(7)
  n <- 60
  covariates <- data.frame(arm = c(rep(1, n / 2), rep(0, n / 2)))
  data <- MCC::GenData(
    covariates = covariates,
    n_event_types = 2L,
    base_event_rate = c(0.8, 1.2)
  )
  expect_error(
    MCC::CompareMvAUCs(data, tau = 3),
    NA
  )
})
