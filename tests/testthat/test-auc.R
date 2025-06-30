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
