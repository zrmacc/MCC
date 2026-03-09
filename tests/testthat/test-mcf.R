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

# -----------------------------------------------------------------------------

test_that("CalcMCF with calc_var = TRUE returns se_mcf and var_mcf.", {
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 1, 0, 1),
    idx = c(1, 2, 3, 4)
  )
  out <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    calc_var = TRUE
  )
  expect_true("se_mcf" %in% names(out))
  expect_true("var_mcf" %in% names(out))
  expect_true(all(out$se_mcf >= 0))
  expect_true(all(out$var_mcf >= 0))
})

test_that("CalcMCF with calc_var = FALSE returns zero variance.", {
  data <- data.frame(
    time = c(1, 2),
    status = c(1, 0),
    idx = c(1, 2)
  )
  out <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    calc_var = FALSE
  )
  # Implementation always returns var_mcf and se_mcf; when calc_var = FALSE they are zero.
  expect_equal(out$var_mcf, c(0, 0))
  expect_equal(out$se_mcf, c(0, 0))
})


# -----------------------------------------------------------------------------

test_that("Test MCF weighting.", {
  
  GetMCF <- function(data){
    CalcMCF(
      idx = data$idx,
      status = data$status,
      time = data$time,
      weights = data$weights
    )
  }
  
  # Case 1.
  data <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 1, 2),
    idx = c(1, 1, 1),
    weights = c(10, 1, 0)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(10, 11, 11))
  
  # Case 2.
  data <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 1, 0),
    idx = c(1, 1, 1),
    weights = c(1, 2, 1)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1, 3, 3))
  
  # Case 3.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 1, 0, 1),
    idx = c(1, 2, 3, 4),
    weights = c(1, 2, 1, 2)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1/4, 3/4, 3/4, 3/4 + 2/3))
  
  # Case 4. 
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 0, 1),
    idx = c(1, 2, 3, 4),
    weights = c(1, 2, 2, 1)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(1/4, 1/4, 1/4, 1/4 + 1/2))
  
  # Case 5.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 2, 1),
    idx = c(1, 2, 3, 4),
    weights = c(2, 1, 0, 2)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(2/4, 2/4, 2/4, 2/4 + 2/2 * 2/3))
  
  # Case 6.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 2, 2, 1),
    idx = c(1, 2, 3, 4),
    weights = c(0, 0, 0, 0)
  )
  observed <- GetMCF(data)
  expect_equal(observed$mcf, c(0, 0, 0, 0))
  
})


# -----------------------------------------------------------------------------
# InfluenceFunction type = "MCF" (MCF influence at tau)
# -----------------------------------------------------------------------------

test_that("InfluenceFunction type MCF returns idx, psi at tau (one row per subject).", {
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 1, 0, 1),
    idx = c(1, 2, 3, 4)
  )
  out <- InfluenceFunction(data, tau = 4, type = "MCF")
  expect_s3_class(out, "data.frame")
  expect_true("idx" %in% names(out))
  expect_true("psi" %in% names(out))
  expect_false("time" %in% names(out))
  expect_equal(nrow(out), length(unique(data$idx)))
})

test_that("InfluenceFunction type MCF variance at tau; var_mcf is monotonicity-adjusted.", {
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 1, 0, 1),
    idx = c(1, 2, 3, 4)
  )
  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    calc_var = TRUE
  )
  psi_df <- InfluenceFunction(data, tau = 4, type = "MCF", mcf = mcf)
  n <- nrow(psi_df)
  var_from_psi <- mean(psi_df$psi^2)
  tau_eff <- max(mcf$time[mcf$time <= 4])
  var_mcf_at_tau <- mcf$var_mcf[mcf$time == tau_eff][1]
  expect_true(var_from_psi <= var_mcf_at_tau + 1e-10)
  psi_first <- InfluenceFunction(data, tau = mcf$time[1], type = "MCF", mcf = mcf)
  expect_equal(mean(psi_first$psi^2), mcf$var_mcf[1])
})

test_that("InfluenceFunction type MCF with mcf = NULL matches precomputed mcf.", {
  withr::local_seed(207)
  data <- GenData(n = 20, tau = 3)
  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    calc_var = FALSE
  )
  out_null <- InfluenceFunction(data, tau = 2, type = "MCF")
  out_mcf <- InfluenceFunction(data, tau = 2, type = "MCF", mcf = mcf)
  expect_equal(out_null$idx, out_mcf$idx)
  expect_equal(out_null$psi, out_mcf$psi)
})


