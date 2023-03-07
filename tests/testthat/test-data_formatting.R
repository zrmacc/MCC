test_that("Expect no change in the case of correctly formatted data.", {
  
  # Correctly formatted data.
  data <- data.frame(
    idx = c(1, 1),
    arm = c(1, 1),
    status = c(1, 0),
    time = c(1, 2),
    strata = c(1, 1)
  )
  
  out <- MCC::FormatData(
    data,
    covars = NULL,
    strata = data$strata,
    cens_after_last = TRUE
  )
  
  expect_equal(data, out)
  
})


# -----------------------------------------------------------------------------


test_that("Check addition of censoring time when missing.", {
  
  # Data without censoring time.
  data <- data.frame(
    idx = c(1, 1),
    arm = c(1, 1),
    status = c(1, 1),
    time = c(1, 2),
    strata = c(1, 1)
  )
  
  # Censor after last set to FALSE.
  observed <- suppressWarnings({
    MCC::FormatData(
      data,
      covars = NULL,
      strata = data$strata,
      cens_after_last = FALSE
    )
  })
  expect_equal(observed, data)
  
  # Censor after last set to TRUE.
  observed <- MCC::FormatData(
      data,
      covars = NULL,
      strata = data$strata,
      cens_after_last = TRUE
  )
  expected <- rbind(data, data[2, ])
  expected$status[3] <- 0
  expect_equal(observed, expected, ignore_attr = TRUE)
})


# -----------------------------------------------------------------------------

test_that("Multiple censoring times triggers an error.", {
  
  data <- data.frame(
    idx = c(1, 1),
    time = c(1, 2),
    status = c(0, 0),
    arm = c(1, 1)
  )
  
  expect_error({
    MCC::FormatData(
      data,
      covars = NULL,
      strata = NULL,
      cens_after_last = TRUE
    )
  })
  
})

# -----------------------------------------------------------------------------

test_that("Character index converted to integer.", {
  
  data <- data.frame(
    idx = c("a", "b"),
    arm = c(1, 1),
    status = c(0, 0),
    time = c(1, 2)
  )
  
  out <- MCC::FormatData(
    data,
    covars = NULL,
    cens_after_last = TRUE
  )
  
  data$idx <- c(1, 2)
  data$strata <- c(1, 1)
  expect_equal(data, out)
  
})


