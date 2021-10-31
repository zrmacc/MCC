# Purpose: Main function for package.
# Updated: 2020-12-12

#' Inference on the Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0).
#' 
#' Two methods of p-value calculation are available. For 'perm', treatment
#' assignments are permuted on each iteration, and the p-value is the
#' proportion of the *null* statistics that are as or more extreme than
#' the *observed* statistics. For 'boot', the p-value is twice the proportion
#' of bootstrap replicates on which the sign of the difference is areas is
#' reversed.
#'
#' @param idx Unique subject index.
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param tau Truncation time.
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
#' @param strata Optional stratification factor.
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @param alpha Alpha level.
#' @param boot Logical, construct bootstrap confidence intervals?
#' @param perm Logical, perform permutation test?
#' @param reps Replicates for bootstrap/permutation inference.
#' @export
#' @return Object of class compAUCs with these slots:
#' \itemize{
#'   \item `@Areas`: The AUC for each arm.
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals.
#'   \item `@Curves`: Mean cumulative count curve for each arm; averaged across strata
#'     if present.
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics.
#'   \item `@Weights`: Per-stratum weights and AUCs.
#' }
#' @examples
#' \donttest{
#' # Simulate data set.
#' covariates <- data.frame(arm = c(rep(1, 50), rep(0, 50)))
#' data <- GenData(
#'   beta_event = log(0.5),
#'   covariates = covariates
#' )
#'
#' aucs <- CompareAUCs(
#'   time = data$time,
#'   status = data$status,
#'   arm = data$arm,
#'   idx = data$idx,
#'   tau = 2,
#'   boot = TRUE,
#'   perm = TRUE,
#'   reps = 25,
#'   alpha = 0.05
#' )
#' show(aucs)
#' }

CompareAUCs <- function(
  idx,
  time,
  status,
  arm,
  tau,
  covars = NULL,
  strata = NULL,
  cens_after_last = TRUE,
  alpha = 0.05,
  boot = FALSE,
  perm = FALSE,
  reps = 2000
) {
  
  if (!is.null(covars) & !is.null(strata)) {
    msg <- paste0(
      "If adjustment for both strata and covariates is needed,\n",
      "include strata indicators within the covariates.")
    stop(msg)
  }
  
  # Format data.
  data <- FormatData(
    idx = idx,
    time = time,
    status = status,
    arm = arm,
    covars = covars,
    strata = strata,
    cens_after_last = cens_after_last
  )
  
  if (!is.null(covars)) {
    out <- CompareAugAUCs(
      data = data,
      tau = tau,
      alpha = alpha,
      boot = boot,
      perm = perm,
      reps = reps
    )
  } else {
    out <- CompareStratAUCs(
      data = data,
      tau = tau,
      alpha = alpha,
      boot = boot,
      perm = perm,
      reps = reps
    )
  }
  
  return(out)
}

