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
  
})
