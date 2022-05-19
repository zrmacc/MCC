test_that("Test tabulation.", {
  
  GetTab <- function(data){
    CalcMCF(
      idx = data$idx,
      status = data$status,
      time = data$time
    )
  }
  
  # Case 1.
  data <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 1, 2),
    idx = c(1, 1, 1)
  )
  observed <- GetTab(data)
  expect_equal(observed$censor, c(0, 0, 0))
  expect_equal(observed$events, c(1, 1, 0))
  expect_equal(observed$death, c(0, 0, 1))
  
  # Case 2.
  data <- data.frame(
    time = c(1, 1, 2),
    status = c(1, 1, 0),
    idx = c(1, 2, 3)
  )
  observed <- GetTab(data)
  expect_equal(observed$censor, c(0, 1))
  expect_equal(observed$events, c(2, 0))
  expect_equal(observed$death, c(0, 0))
  
  # Case 3.
  data <- data.frame(
    time = c(1, 1, 2),
    status = c(2, 2, 0),
    idx = c(1, 2, 3)
  )
  observed <- GetTab(data)
  expect_equal(observed$censor, c(0, 1))
  expect_equal(observed$events, c(0, 0))
  expect_equal(observed$death, c(2, 0))
  expect_equal(observed$nar, c(3, 1))
  
  # Case 4.
  data <- data.frame(
    time = c(1, 1, 2, 2, 3, 3),
    status = c(1, 2, 1, 0, 1, 1),
    idx = c(1, 1, 2, 2, 3, 4)
  )
  observed <- GetTab(data)
  expect_equal(observed$censor, c(0, 1, 0))
  expect_equal(observed$events, c(1, 1, 2))
  expect_equal(observed$death, c(1, 0, 0))
  expect_equal(observed$nar, c(4, 3, 2))

})

# -----------------------------------------------------------------------------

test_that("Test MCF calculation.", {
  
  
  GetMCF <- function(data){
    CalcMCF(
      idx = data$idx,
      status = data$status,
      time = data$time
    )
  }
  
  # Case 1.
  data <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 1, 2),
    idx = c(1, 1, 1)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1, 2, 2))
  
  # Case 2.
  data <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 1, 0),
    idx = c(1, 1, 1)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1, 2, 2))
  
  # Case 3.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 1, 0, 1),
    idx = c(1, 2, 3, 4)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1/4, 2/4, 2/4, 2/4 + 1/3))
  
  # Case 4. 
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 0, 1),
    idx = c(1, 2, 3, 4)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1/4, 1/4, 1/4, 1/4 + 1/2))
  
  # Case 5.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 2, 1),
    idx = c(1, 2, 3, 4)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1/4, 1/4, 1/4, 1/4 + 2/3 * 1/2))
  
  # Case 6.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 2, 2, 1),
    idx = c(1, 2, 3, 4)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1/4, 1/4, 1/4, 1/4 + 1/2 * 1/2))
  
})