# Tests for process-specific eligibility influence-function weighting (Derivations §1.5).

#' Build simple multivariate long-format data for testing.
#'
#' @param specs List of per-subject specifications.
#' @return Long-format data.frame.
#' @noRd
.MvTestData <- function(specs) {
  rows <- lapply(specs, function(s) {
    arm <- s$arm
    idx <- s$idx
    out <- NULL
    for (ev in s$events) {
      out <- rbind(
        out,
        data.frame(
          arm = arm,
          idx = idx,
          time = ev$time,
          status = ev$status,
          event_type = ev$event_type,
          stringsAsFactors = FALSE
        )
      )
    }
    return(out)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}


#' Partial-eligibility toy data with two-thirds eligible per process.
#'
#' @return Long-format multivariate data.frame.
.PartialEligTestData <- function() {
  .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 1, idx = 2,
      events = list(
        list(time = 3, status = 0, event_type = 1L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 1, idx = 3,
      events = list(
        list(time = 1, status = 1, event_type = 2L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 4,
      events = list(
        list(time = 2, status = 1, event_type = 1L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 5,
      events = list(
        list(time = 3, status = 0, event_type = 1L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 6,
      events = list(
        list(time = 1, status = 1, event_type = 2L),
        list(time = 3, status = 0, event_type = NA)
      )
    )
  ))
}

#' Recompute arm-theta components without eligibility scaling (old bug).
#'
#' @param arm_data Arm-specific formatted data.
#' @param event_types Event types.
#' @param tau Truncation time.
#' @return List with psi_mat and n_arm.
.UnweightedMvArmPsi <- function(arm_data, event_types, tau) {
  all_idx <- sort(unique(arm_data$idx))
  n_arm <- length(all_idx)
  K <- length(event_types)
  psi_mat <- matrix(0, nrow = n_arm, ncol = K)
  rownames(psi_mat) <- as.character(all_idx)
  colnames(psi_mat) <- as.character(event_types)
  for (k in seq_len(K)) {
    et <- event_types[k]
    at_risk_idx <- MCC:::.MvAtRiskIdx(data = arm_data, k = et)
    arm_k <- arm_data[arm_data$idx %in% at_risk_idx, , drop = FALSE]
    proc <- MCC:::.SubsetProcessK(
      data = arm_k,
      k = et,
      cens_after_last = TRUE
    )
    fit <- MCC:::MvProcessAUC(data = proc, tau = tau, calc_var = TRUE)
    row_match <- match(as.character(fit$psi$idx), rownames(psi_mat))
    psi_mat[row_match, k] <- fit$psi$psi
  }
  out <- list(psi_mat = psi_mat, n_arm = n_arm)
  return(out)
}


test_that("Eligible subjects receive Q/pi weight on influence functions.", {
  data <- .PartialEligTestData()
  fmt <- MCC:::FormatMvData(data = data, event_type_name = "event_type")
  arm1 <- fmt$data[fmt$data$arm == 1, , drop = FALSE]
  theta <- MCC:::MvArmTheta(
    arm_data = arm1,
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL
  )
  at_risk_k1 <- MCC:::.MvAtRiskIdx(data = arm1, k = 1L)
  arm_k1 <- arm1[arm1$idx %in% at_risk_k1, , drop = FALSE]
  proc_k1 <- MCC:::.SubsetProcessK(
    data = arm_k1,
    k = 1L,
    cens_after_last = TRUE
  )
  raw_k1 <- MCC:::MvProcessAUC(data = proc_k1, tau = 3, calc_var = TRUE)
  pi_hat <- raw_k1$n / theta$n_arm
  weight <- 1 / pi_hat
  eligible_rows <- as.character(raw_k1$psi$idx)
  ineligible_rows <- setdiff(rownames(theta$psi_mat), eligible_rows)
  expect_equal(
    unname(theta$psi_mat[eligible_rows, "1"]),
    unname(weight * raw_k1$psi$psi),
    tolerance = 1e-10
  )
  expect_true(all(theta$psi_mat[ineligible_rows, "1"] == 0))
  expect_equal(theta$psi_mat[ineligible_rows, "2"], 0)
})


test_that("Arm covariance matches (1/n_a) sum phi phi' with eligibility weights.", {
  data <- .PartialEligTestData()
  fmt <- MCC:::FormatMvData(data = data, event_type_name = "event_type")
  arm0 <- fmt$data[fmt$data$arm == 0, , drop = FALSE]
  theta0 <- MCC:::MvArmTheta(
    arm_data = arm0,
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL
  )
  sigma_manual <- crossprod(theta0$psi_mat) / theta0$n_arm
  obs <- MCC:::CalcMvAugAUC(
    data = fmt$data,
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL,
    return_areas = FALSE
  )
  expect_equal(obs$cov_theta[["0"]], sigma_manual, tolerance = 1e-10)
})


test_that("Contrast variance matches Sigma_1/n_1 + Sigma_0/n_0.", {
  data <- .PartialEligTestData()
  fmt <- MCC:::FormatMvData(data = data, event_type_name = "event_type")
  obs <- MCC:::CalcMvAugAUC(
    data = fmt$data,
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL,
    return_areas = FALSE
  )
  theta0 <- MCC:::MvArmTheta(
    arm_data = fmt$data[fmt$data$arm == 0, , drop = FALSE],
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL
  )
  theta1 <- MCC:::MvArmTheta(
    arm_data = fmt$data[fmt$data$arm == 1, , drop = FALSE],
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL
  )
  sigma0 <- crossprod(theta0$psi_mat) / theta0$n_arm
  sigma1 <- crossprod(theta1$psi_mat) / theta1$n_arm
  cov_manual <- sigma1 / theta1$n_arm + sigma0 / theta0$n_arm
  expect_equal(obs$cov_delta, cov_manual, tolerance = 1e-10)
})


test_that("Full eligibility leaves influence weights at unity.", {
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 1, event_type = 2L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 1, event_type = 2L),
        list(time = 3, status = 0, event_type = NA)
      )
    )
  ))
  fmt <- MCC:::FormatMvData(data = data, event_type_name = "event_type")
  arm1 <- fmt$data[fmt$data$arm == 1, , drop = FALSE]
  theta <- MCC:::MvArmTheta(
    arm_data = arm1,
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL
  )
  unweighted <- .UnweightedMvArmPsi(
    arm_data = arm1,
    event_types = c(1L, 2L),
    tau = 3
  )
  expect_equal(theta$psi_mat, unweighted$psi_mat, tolerance = 1e-10)
})


test_that("Omitting eligibility weight underestimates contrast SE.", {
  data <- .PartialEligTestData()
  fmt <- MCC:::FormatMvData(data = data, event_type_name = "event_type")
  obs <- MCC:::CalcMvAugAUC(
    data = fmt$data,
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = NULL,
    return_areas = FALSE
  )
  theta0 <- .UnweightedMvArmPsi(
    arm_data = fmt$data[fmt$data$arm == 0, , drop = FALSE],
    event_types = c(1L, 2L),
    tau = 3
  )
  theta1 <- .UnweightedMvArmPsi(
    arm_data = fmt$data[fmt$data$arm == 1, , drop = FALSE],
    event_types = c(1L, 2L),
    tau = 3
  )
  sigma0 <- crossprod(theta0$psi_mat) / theta0$n_arm
  sigma1 <- crossprod(theta1$psi_mat) / theta1$n_arm
  cov_buggy <- sigma1 / theta1$n_arm + sigma0 / theta0$n_arm
  se_buggy <- sqrt(diag(cov_buggy)[1])
  se_correct <- sqrt(diag(obs$cov_delta)[1])
  expect_lt(se_buggy, se_correct)
})


test_that("Augmentation uses process-specific eligible covariate means.", {
  data <- .PartialEligTestData()
  data$X <- c(0, 1, 2, -1, 0.5, 1.5)
  fmt <- MCC:::FormatMvData(
    data = data,
    event_type_name = "event_type",
    covars = data[, "X", drop = FALSE]
  )
  theta1 <- MCC:::MvArmTheta(
    arm_data = fmt$data[fmt$data$arm == 1, , drop = FALSE],
    event_types = c(1L, 2L),
    tau = 3,
    cens_after_last = TRUE,
    covar_cols = "X"
  )
  arm1 <- fmt$data[fmt$data$arm == 1, , drop = FALSE]
  at_risk_k1 <- MCC:::.MvAtRiskIdx(data = arm1, k = 1L)
  at_risk_k2 <- MCC:::.MvAtRiskIdx(data = arm1, k = 2L)
  xbar_k1 <- mean(arm1$X[arm1$idx %in% at_risk_k1])
  xbar_k2 <- mean(arm1$X[arm1$idx %in% at_risk_k2])
  expect_equal(unname(theta1$xbar_by_k[[1]]), unname(xbar_k1), tolerance = 1e-10)
  expect_equal(unname(theta1$xbar_by_k[[2]]), unname(xbar_k2), tolerance = 1e-10)
})


test_that("Thirds eligibility DGP yields calibrated Wald SEs.", {
  set.seed(20260623)
  reps <- 250L
  deltas <- numeric(reps)
  ses <- numeric(reps)
  for (i in seq_len(reps)) {
    data <- SimMvEligData(
      scenario_name = "elig_no_covar",
      study = "t1e",
      n_per_arm = 50L,
      tau_gen = 4
    )
    mv <- CompareMvAUCs(data = data, tau = 2, reps = NULL, covars = NULL)
    cis <- mv@CIs[mv@CIs$method == "asymptotic" & mv@CIs$event_type == 1L, ]
    deltas[i] <- cis$observed
    ses[i] <- cis$se
  }
  ratio <- sd(deltas) / mean(ses)
  expect_gt(ratio, 0.80)
  expect_lt(ratio, 1.20)
})
