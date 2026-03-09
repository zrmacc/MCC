# InfluenceFunction (idx, psi at tau for MCF or AUC)
# -----------------------------------------------------------------------------

test_that("InfluenceFunction returns data.frame with idx, psi (one row per subject).", {
  data <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3),
    time = c(1, 2, 1.5, 2.5, 0.5, 2),
    status = c(1, 0, 1, 0, 1, 0)
  )
  for (typ in c("AUC", "MCF")) {
    out <- InfluenceFunction(data, tau = 2, type = typ)
    expect_s3_class(out, "data.frame")
    expect_true("idx" %in% names(out))
    expect_true("psi" %in% names(out))
    expect_equal(nrow(out), length(unique(data$idx)))
  }
})

test_that("InfluenceFunction type = 'AUC' matches PsiAUC result.", {
  withr::local_seed(208)
  data <- GenData(n = 25, tau = 3)
  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    calc_var = FALSE
  )
  out_inf <- InfluenceFunction(data, tau = 2, type = "AUC", mcf = mcf)
  psi_auc <- MCC:::PsiAUC(
    event_rate = mcf$weighted_event_rate,
    grid_time = mcf$time,
    idx = data$idx,
    haz = mcf$haz,
    nar = mcf$nar,
    status = data$status,
    surv = mcf$surv,
    tau = 2,
    time = data$time,
    weights = rep(1, nrow(data))
  )
  expect_equal(out_inf$idx, psi_auc$idx)
  expect_equal(out_inf$psi, psi_auc$psi)
})

test_that("InfluenceFunction type = 'MCF' returns idx, psi consistent with MCF at tau.", {
  data <- data.frame(
    idx = c(1, 2, 3),
    time = c(1, 1, 2),
    status = c(1, 1, 0)
  )
  mcf <- CalcMCF(idx = data$idx, status = data$status, time = data$time, calc_var = FALSE)
  out <- InfluenceFunction(data, tau = 2, type = "MCF", mcf = mcf)
  expect_equal(nrow(out), 3)
  tau_eff <- max(mcf$time[mcf$time <= 2])
  mcf_at_tau <- mcf$mcf[mcf$time == tau_eff][1]
  pseudo <- mcf_at_tau + out$psi
  expect_true(all(is.finite(pseudo)))
})

test_that("InfluenceFunction with mcf = NULL matches precomputed mcf.", {
  data <- data.frame(
    idx = c(1, 2, 3),
    time = c(1, 1, 2),
    status = c(1, 1, 0)
  )
  out_null <- InfluenceFunction(data, tau = 2, type = "MCF")
  mcf <- CalcMCF(idx = data$idx, status = data$status, time = data$time, calc_var = FALSE)
  out_mcf <- InfluenceFunction(data, tau = 2, type = "MCF", mcf = mcf)
  expect_equal(out_null$idx, out_mcf$idx)
  expect_equal(out_null$psi, out_mcf$psi)
})

test_that("InfluenceFunction invalid type errors.", {
  data <- data.frame(idx = 1, time = 1, status = 0)
  expect_error(InfluenceFunction(data, tau = 1, type = "invalid"), "arg")
})
