test_that("Reused subject IDs raises error.", {
  data <- data.frame(
    time = c(1, 1),
    status = c(1, 2),
    idx = c(1, 1),
    arm = c(1, 0),
    weights = c(1, 1)
  )
  
  expect_error(CompareAUCs(data, tau = 1))
})

# -----------------------------------------------------------------------------
# Base.
# -----------------------------------------------------------------------------

test_that("Test AUC comparison.", {
  
  # Truncation time.
  tau <- 4
  
  GetAUCDiff <- function(data){
    out <-suppressWarnings(
      CompareAUCs(
        data,
        tau = tau,
        cens_after_last = FALSE
      ))
    return(out)
  }
  
  # Case 1.
  data1 <- data.frame(
    time = c(1, 1, 2, 2),
    status = c(1, 1, 1, 1),
    idx = c(1, 2, 3, 4)
  )
  data1$arm <- 1
  auc1 <- 0.5 * (2 - 1) + 1.0 * (tau - 2)
  
  data0 <- data.frame(
    time = c(1, 1, 2, 2),
    status = c(0, 2, 0, 2),
    idx = c(5, 6, 7, 8)
  )
  data0$arm <- 0
  auc0 <- 0
  
  data <- rbind(data1, data0)
  
  observed <- GetAUCDiff(data)
  expect_equal(observed@MargAreas$area[1], auc0)
  expect_equal(observed@MargAreas$area[2], auc1)
  expect_equal(observed@Pvals$observed[1], auc1 - auc0)
  
  
  # Case 2.
  data1 <- data.frame(
    time = c(1, 2, 3, 4),
    status =c(1, 2, 2, 1),
    idx = c(1, 2, 3, 4)
  )
  data1$arm <- 1
  auc1 <- 0.25 * (4 - 1) + 0.5 * (tau - 4)
  
  data0 <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(2, 2, 1, 0),
    idx = c(5, 6, 7, 8)
  )
  data0$arm <- 0
  auc0 <- 0.25 * (tau - 3)
  
  data <- rbind(data1, data0)
  
  observed <- GetAUCDiff(data)
  expect_equal(observed@MargAreas$area[1], auc0)
  expect_equal(observed@MargAreas$area[2], auc1)
  expect_equal(observed@Pvals$observed[1], auc1 - auc0)

})


# -----------------------------------------------------------------------------
# Augmentation.
# -----------------------------------------------------------------------------

test_that("Overall test of augmentation analysis.", {
  n <- 100
  covariates <- data.frame(
    arm = c(rep(1, n/2), rep(0, n/2)),
    covar = rnorm(n)
  )
  data <- GenData(
    beta_event = c(log(0.5), 1),
    covariates = covariates
  )

  expect_error(CompareAUCs(
    data,
    covars = data$covar,
    tau = 2,
    boot = FALSE,
    perm = FALSE,
    reps = 25,
    alpha = 0.05
  ), NA)
})

# -----------------------------------------------------------------------------
# Augmentation.
# -----------------------------------------------------------------------------

test_that("Overall test of augmentation analysis.", {
  n <- 100
  covariates <- data.frame(
    arm = c(rep(1, n/2), rep(0, n/2)),
    covar = rnorm(n)
  )
  data <- GenData(
    beta_event = c(log(0.5), 1),
    covariates = covariates
  )
  
  Run <- function(){
    CompareAUCs(
      data,
      covars = data$covar,
      tau = 2,
      boot = TRUE,
      perm = TRUE,
      reps = 25,
      alpha = 0.05
    )
  }
  
  expect_error(Run(), NA)
})

# -----------------------------------------------------------------------------
# Stratification.
# -----------------------------------------------------------------------------

test_that("Overall test of stratified analysis.", {
  n <- 100
  covariates <- data.frame(
    arm = c(rep(1, n/2), rep(0, n/2)),
    strata = stats::rbinom(n, 1, 0.2)
  )
  data <- GenData(
    beta_event = c(log(0.5), log(0.8)),
    covariates = covariates
  )
  
  Run <- function() {
    CompareAUCs(
      data,
      strata = data$strata,
      tau = 2,
      boot = TRUE,
      perm = TRUE,
      reps = 25,
      alpha = 0.05
    )
  }
  
  # Expect no error.
  expect_error(Run(), NA)
})

# -----------------------------------------------------------------------------
# SingleArmAUC
# -----------------------------------------------------------------------------

test_that("SingleArmAUC returns object with expected slots and area matches AUC.", {
  # Each subject has a terminal event (censor or death) so no warning.
  data <- data.frame(
    time = c(1, 2, 2, 2, 3, 3, 3, 4, 4),
    status = c(1, 1, 0, 1, 0, 1, 0, 1, 0),
    idx = c(1, 1, 1, 2, 2, 3, 3, 4, 4)
  )
  tau <- 4
  out <- SingleArmAUC(data, tau = tau, boot = FALSE, cens_after_last = FALSE)
  expect_s4_class(out, "CompStratAUCs")
  expect_true("MargAreas" %in% slotNames(out))
  area_from_slot <- out@MargAreas$area[1]
  mcf <- CalcMCF(idx = data$idx, status = data$status, time = data$time, calc_var = FALSE)
  area_from_auc <- AUC(mcf$time, mcf$mcf, tau)
  expect_equal(area_from_slot, area_from_auc, tolerance = 1e-9)
})

# -----------------------------------------------------------------------------
# CalcMargMCF
# -----------------------------------------------------------------------------

test_that("CalcMargMCF returns data.frame with arm, time, mcf for two-arm stratified data.", {
  withr::local_seed(301)
  covariates <- data.frame(
    arm = c(rep(1, 40), rep(0, 40)),
    strata = rep(c(0, 1), 40)
  )
  data <- GenData(
    beta_event = c(log(0.8), log(1.2)),
    covariates = covariates,
    tau = 3
  )
  data$weights <- 1
  out <- CalcMargMCF(data)
  expect_s3_class(out, "data.frame")
  expect_true("arm" %in% names(out))
  expect_true("time" %in% names(out))
  expect_true("mcf" %in% names(out))
  expect_equal(sort(unique(out$arm)), c(0, 1))
})

