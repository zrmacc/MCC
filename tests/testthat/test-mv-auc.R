# Tests for CompareMvAUCs and multivariate AUMCF estimation.

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

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

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  return(x)
}


# -----------------------------------------------------------------------------
# K = 1 equivalence
# -----------------------------------------------------------------------------

test_that("K = 1 matches CompareAUCs without covariates.", {
  tau <- 4

  data1 <- data.frame(
    time = c(1, 1, 1, 1, 2, 2, 2, 2),
    status = c(1, 0, 1, 0, 1, 0, 1, 0),
    idx = rep(1:4, each = 2),
    arm = 1
  )
  data0 <- data.frame(
    time = c(1, 1, 2, 2),
    status = c(0, 2, 0, 2),
    idx = c(5, 6, 7, 8),
    arm = 0
  )
  data <- rbind(data1, data0)

  strat <- CompareAUCs(
    data,
    tau = tau,
    cens_after_last = FALSE
  )
  mv <- CompareMvAUCs(
    data,
    tau = tau,
    cens_after_last = FALSE
  )

  expect_equal(mv@Areas$area, strat@MargAreas$area, tolerance = 1e-10)
  expect_equal(
    mv@Pvals$observed[mv@Pvals$event_type == 1],
    strat@Pvals$observed[strat@Pvals$contrast == "A1-A0"],
    tolerance = 1e-10
  )
  se_strat <- sqrt(sum(strat@MargAreas$se^2))
  expect_equal(sqrt(mv@CovDelta[1, 1]), se_strat, tolerance = 1e-8)
})


test_that("K = 1 matches CompareAUCs with covariates.", {
  set.seed(42)
  n <- 100
  covariates <- data.frame(
    arm = c(rep(1, n / 2), rep(0, n / 2)),
    covar = rnorm(n)
  )
  data <- GenData(
    beta_event = c(log(0.5), 1),
    covariates = covariates
  )

  aug <- CompareAUCs(
    data,
    covars = data$covar,
    tau = 2
  )
  mv <- CompareMvAUCs(
    data,
    covars = data$covar,
    tau = 2
  )

  expect_equal(mv@Areas$area, aug@Areas$area, tolerance = 1e-10)
  expect_equal(
    mv@Pvals$observed[mv@Pvals$event_type == 1],
    aug@Pvals$observed[aug@Pvals$contrast == "A1-A0"],
    tolerance = 1e-10
  )
  expect_equal(
    sqrt(mv@CovDelta[1, 1]),
    aug@CIs$se[aug@CIs$contrast == "A1-A0" & aug@CIs$method == "asymptotic"],
    tolerance = 1e-8
  )
})


test_that("K = 1 bootstrap matches asymptotic observed statistics.", {
  set.seed(7)
  n <- 80
  covariates <- data.frame(
    arm = c(rep(1, n / 2), rep(0, n / 2)),
    covar = rnorm(n)
  )
  data <- GenData(
    beta_event = c(log(0.5), 1),
    covariates = covariates
  )

  mv_asym <- CompareMvAUCs(data, covars = data$covar, tau = 2)
  mv_boot <- CompareMvAUCs(
    data,
    covars = data$covar,
    tau = 2,
    reps = 30
  )

  expect_equal(
    mv_boot@Pvals$observed[mv_boot@Pvals$method == "asymptotic"],
    mv_asym@Pvals$observed,
    tolerance = 1e-10
  )
  expect_true(any(mv_boot@Pvals$method == "bootstrap"))
})


# -----------------------------------------------------------------------------
# Multivariate deterministic
# -----------------------------------------------------------------------------

test_that("Multivariate areas match hand calculations for K = 2.", {
  tau <- 4
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
        list(time = 2, status = 2, event_type = 2L),
        list(time = 2, status = 2, event_type = NA)
      )
    )
  ))

  mv <- CompareMvAUCs(data, tau = tau, cens_after_last = FALSE)

  expect_equal(
    mv@Areas$area[mv@Areas$arm == 1 & mv@Areas$event_type == 1],
    3
  )
  expect_equal(
    mv@Areas$area[mv@Areas$arm == 0 & mv@Areas$event_type == 1],
    3
  )
  expect_equal(
    mv@Pvals$observed[mv@Pvals$event_type == 1],
    0,
    tolerance = 1e-10
  )
  expect_equal(
    mv@Pvals$observed[mv@Pvals$event_type == 2],
    2,
    tolerance = 1e-10
  )
  expect_equal(length(unique(mv@Pvals$event_type)), 2)
})


# -----------------------------------------------------------------------------
# Risk sets
# -----------------------------------------------------------------------------

test_that("Risk set inferred from event_type rows changes sample size.", {
  tau <- 3
  base <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 1, event_type = 2L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 1, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 3,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 4,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  mv <- CompareMvAUCs(
    base,
    tau = tau,
    cens_after_last = FALSE
  )

  expect_equal(
    mv@Areas$n[mv@Areas$arm == 1 & mv@Areas$event_type == 2],
    1
  )
  expect_equal(
    mv@Areas$n[mv@Areas$arm == 0 & mv@Areas$event_type == 2],
    1
  )
})


test_that("Typed terminal only places subject at risk with zero events.", {
  tau <- 4
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = 1L),
        list(time = 3, status = 0, event_type = 2L)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 0, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  mv <- CompareMvAUCs(data, tau = tau, cens_after_last = FALSE)
  area1_t2 <- mv@Areas$area[mv@Areas$arm == 1 & mv@Areas$event_type == 2]
  n1_t2 <- mv@Areas$n[mv@Areas$arm == 1 & mv@Areas$event_type == 2]
  expect_equal(n1_t2, 1)
  expect_equal(area1_t2, 0)
})


# -----------------------------------------------------------------------------
# Typed vs global terminal events
# -----------------------------------------------------------------------------

test_that("Typed terminal ends one process only.", {
  tau <- 5
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = 1L),
        list(time = 3, status = 1, event_type = 2L),
        list(time = 5, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  mv <- CompareMvAUCs(data, tau = tau, cens_after_last = FALSE)
  area1_t2 <- mv@Areas$area[mv@Areas$arm == 1 & mv@Areas$event_type == 2]
  expect_equal(area1_t2, 2)
})


test_that("Partial typed terminals receive cens_after_last for other types.", {
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = 1L),
        list(time = 3, status = 1, event_type = 2L)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  proc <- MCC:::.SubsetProcessK(
    data = data,
    k = 2L,
    cens_after_last = TRUE
  )
  subj1 <- proc[proc$idx == 1, , drop = FALSE]
  expect_equal(subj1$time, c(3, 3))
  expect_equal(subj1$status, c(1L, 0L))

  expect_warning(
    MCC:::FormatMvData(
      data = data,
      cens_after_last = FALSE,
      event_type_name = "event_type"
    ),
    "Patients without observation terminating events were found for: idx 1, event_type 2"
  )
})


test_that("Global terminal ends follow-up for all event types.", {
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 1, event_type = 2L),
        list(time = 4, status = 0, event_type = NA)
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

  proc1 <- MCC:::.SubsetProcessK(data, k = 1L, cens_after_last = TRUE)
  proc2 <- MCC:::.SubsetProcessK(data, k = 2L, cens_after_last = TRUE)
  expect_equal(
    max(proc1$time[proc1$idx == 1 & proc1$status %in% c(0L, 2L)]),
    4
  )
  expect_equal(
    max(proc2$time[proc2$idx == 1 & proc2$status %in% c(0L, 2L)]),
    4
  )
  expect_no_warning(
    CompareMvAUCs(data, tau = 5, cens_after_last = FALSE)
  )
})


test_that("Partial eligibility scales marginal influence functions.", {
  data <- .MvTestData(list(
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

  mv <- CompareMvAUCs(data, tau = 3)
  proc <- do.call(rbind, lapply(c(0L, 1L), function(a) {
    arm_data <- data[data$arm == a, , drop = FALSE]
    at_risk_idx <- MCC:::.MvAtRiskIdx(data = arm_data, k = 1L)
    arm_data <- arm_data[arm_data$idx %in% at_risk_idx, , drop = FALSE]
    out <- MCC:::.SubsetProcessK(
      data = arm_data,
      k = 1L,
      cens_after_last = TRUE
    )
    out$arm <- a
    return(out)
  }))
  strat <- CompareAUCs(proc, tau = 3)
  strat_se <- strat@CIs$se[
    strat@CIs$method == "asymptotic" &
      strat@CIs$contrast == "A1-A0"
  ]

  expect_equal(mv@Areas$n[mv@Areas$event_type == 1L], c(2, 2))
  expect_equal(sqrt(mv@CovDelta[1, 1]), strat_se, tolerance = 1e-8)
})


test_that("Multiple typed terminals for one event type trigger an error.", {
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = 1L),
        list(time = 3, status = 2, event_type = 1L)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  expect_error(
    MCC:::.SubsetProcessK(data, k = 1L, cens_after_last = TRUE),
    "multiple observation terminating events for event type 1"
  )
})


test_that("Subjects without type-k rows skip terminal validation for k.", {
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L)
      )
    ),
    list(
      arm = 1, idx = 3,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  expect_warning(
    CompareMvAUCs(
      data,
      tau = 4,
      cens_after_last = FALSE
    ),
    "Patients without observation terminating events were found for: idx 1, event_type 1"
  )
  expect_no_warning(
    CompareMvAUCs(
      data,
      tau = 4,
      cens_after_last = TRUE
    )
  )
})


test_that("Subjects not at risk for event type skip terminal validation.", {
  data <- .MvTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = 1L)
      )
    ),
    list(
      arm = 1, idx = 3,
      events = list(
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  expect_no_warning(
    CompareMvAUCs(
      data,
      tau = 5,
      cens_after_last = FALSE
    )
  )
})


# -----------------------------------------------------------------------------
# Covariance and weights
# -----------------------------------------------------------------------------

test_that("Covariance matrices are symmetric with positive diagonal.", {
  set.seed(11)
  n <- 60
  covariates <- data.frame(
    arm = c(rep(1, n / 2), rep(0, n / 2)),
    covar = rnorm(n)
  )
  data <- GenData(covariates = covariates)
  data$event_type <- 1L

  mv <- CompareMvAUCs(data, tau = 2, cens_after_last = TRUE)
  expect_true(isSymmetric(mv@CovDelta))
  expect_true(isSymmetric(mv@CovTheta[["0"]]))
  expect_true(all(diag(mv@CovDelta) > 0))
})


test_that("Optional weights produce scalar weighted contrast.", {
  tau <- 4
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
        list(time = 1, status = 1, event_type = 2L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))

  mv <- CompareMvAUCs(
    data,
    tau = tau,
    cens_after_last = FALSE,
    process_weights = c(1, 1)
  )

  comp_sum <- sum(mv@Pvals$observed[mv@Pvals$method == "asymptotic"])
  expect_equal(nrow(mv@Weighted), 1)
  expect_equal(mv@Weighted$observed, comp_sum, tolerance = 1e-10)
  expect_equal(nrow(mv@CIs), 2)
  expect_equal(nrow(mv@Pvals), 2)
})


test_that("Augmentation analysis runs and adjusts estimates.", {
  set.seed(99)
  n <- 120
  covariates <- data.frame(
    arm = c(rep(1, n / 2), rep(0, n / 2)),
    covar = c(rnorm(n / 2, 1), rnorm(n / 2, 0))
  )
  data <- GenData(
    beta_event = c(log(0.5), 0.5),
    covariates = covariates
  )
  data$event_type <- 1L

  unadj <- CompareMvAUCs(data, tau = 2)
  adj <- CompareMvAUCs(data, tau = 2, covars = data$covar)

  expect_false(isTRUE(all.equal(
    unadj@Pvals$observed,
    adj@Pvals$observed
  )))
})
