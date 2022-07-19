test_that("Test sampling-based p-values.", {
  
  rej <- c(rep(1, 55), rep(0, 44))
  obs <- CalcP(rej)
  expect_equal(obs, 2 * 45 / 100)
  
  rej <- c(rep(1, 39), rep(0, 60))
  obs <- CalcP(rej)
  expect_equal(obs, 2 * 40 / 100)
  
  rej <- c(rep(1, 49), rep(0, 50))
  obs <- CalcP(rej)
  expect_equal(obs, 2 * 50 / 100)
  
  rej <- rep(0, 99)
  obs <- CalcP(rej)
  expect_equal(obs, 2 * 1 / 100)  
  
})
