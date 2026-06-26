# Tests for multivariate eligibility simulation DGP utilities.

test_that("SplitThirds sums to n.", {
  expect_equal(sum(SplitThirds(50)), 50)
  expect_equal(SplitThirds(50), c(17L, 17L, 16L))
})

test_that("AssignMvEligibility assigns thirds within each arm.", {
  elig <- AssignMvEligibility(n_per_arm = 60)
  for (a in c(0L, 1L)) {
    arm_elig <- elig[elig$arm == a, , drop = FALSE]
    expect_equal(nrow(arm_elig), 60)
    expect_equal(unname(as.vector(table(arm_elig$elig_pattern))), rep(20L, 3L))
  }
})

test_that("ApplyMvEventEligibility enforces process-specific risk sets.", {
  set.seed(1)
  data <- SimMvEligData(
    scenario_name = "elig_no_covar",
    study = "t1e",
    n_per_arm = 30,
    tau_gen = 4
  )
  elig_map <- unique(data[, c("idx", "elig_pattern")])
  k1_only <- elig_map$idx[elig_map$elig_pattern == "k1_only"]
  k2_only <- elig_map$idx[elig_map$elig_pattern == "k2_only"]
  expect_false(any(data$status == 1 & data$idx %in% k1_only & data$event_type == 2))
  expect_false(any(data$status == 1 & data$idx %in% k2_only & data$event_type == 1))
  mv <- CompareMvAUCs(data = data, tau = 3, reps = NULL, covars = NULL)
  n1 <- mv@Areas$n[mv@Areas$event_type == 1]
  n2 <- mv@Areas$n[mv@Areas$event_type == 2]
  expect_true(all(n1 >= 18))
  expect_true(all(n2 >= 18))
})

test_that("CompareMvAUCs risk-set counts match eligibility assignment.", {
  set.seed(2)
  data <- SimMvEligData(
    scenario_name = "elig_no_covar",
    study = "t1e",
    n_per_arm = 50,
    tau_gen = 4
  )
  elig_map <- unique(data[, c("idx", "arm", "elig_pattern")])
  mv <- CompareMvAUCs(data = data, tau = 2, reps = NULL, covars = NULL)
  for (a in c(0L, 1L)) {
    arm_map <- elig_map[elig_map$arm == a, , drop = FALSE]
    actual_k1 <- sum(arm_map$elig_pattern %in% c("both", "k1_only"))
    actual_k2 <- sum(arm_map$elig_pattern %in% c("both", "k2_only"))
    inferred_k1 <- mv@Areas$n[mv@Areas$arm == a & mv@Areas$event_type == 1L]
    inferred_k2 <- mv@Areas$n[mv@Areas$arm == a & mv@Areas$event_type == 2L]
    expect_equal(inferred_k1, actual_k1)
    expect_equal(inferred_k2, actual_k2)
  }
})

test_that("elig_covar_types encodes type-specific beta_event matrix.", {
  cfg <- .MvEligScenarioConfig("elig_covar_types", study = "t1e")
  beta <- cfg$beta_event
  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(4L, 2L))
  expect_equal(beta["X1", "k1"], log(1.25))
  expect_equal(beta["X1", "k2"], 0)
  expect_equal(beta["X2", "k1"], 0)
  expect_equal(beta["X2", "k2"], log(0.75))
})
