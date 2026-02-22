# Purpose: Main function for package.
# Updated: 2024-02-19

# -----------------------------------------------------------------------------
# Single-sample.
# -----------------------------------------------------------------------------

#' Single Arm Area Under the Cumulative Count Curve
#'
#' @param data Data.frame. 
#' @param alpha Type I error level.
#' @param boot Logical, construct bootstrap confidence intervals (slow)?
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @param idx_name Name of column containing a unique subject index.
#' @param reps Number of replicates for bootstrap inference.
#' @param status_name Name of column containing the status. Must be coded as 0 for
#'   censoring, 1 for event, 2 for death. Each subject should have an 
#'   observation-terminating event, either censoring or death.
#' @param strata Optional stratification factor.
#' @param tau Numeric truncation time.
#' @param time_name Name of column containing the observation time.
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Object of class \code{CompStratAUCs} with these slots:
#' \itemize{
#'   \item `@Areas`: The AUC for each arm.
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals.
#'   \item `@Curves`: Mean cumulative count curve for each arm; averaged across strata
#'     if present.
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics.
#'   \item `@Weights`: Per-stratum weights and AUCs.
#' }
#' @export
#' @examples
#' # Simulate data set.
#' covar <- data.frame(strata = rep(c(1, 2), each = 50))
#' data <- GenData(beta_event = log(0.5), covariates = covar)
#' # Calculate AUC.
#' auc <- SingleArmAUC(data, strata = data$strata, tau = 2)
#' \donttest{auc <- SingleArmAUC(data, boot = TRUE, reps = 100, strata = data$strata, tau = 2)}
#' show(auc)
SingleArmAUC <- function(
    data,
    alpha = 0.05,
    boot = FALSE,
    cens_after_last = TRUE,
    idx_name = "idx",
    reps = 2000,
    status_name = "status",
    strata = NULL,
    tau = NULL,
    time_name = "time",
    weights = NULL
) {
  
  # Rename columns as necessary.
  idx <- time <- status <- NULL
  data <- data %>%
    dplyr::rename(
      idx = {{idx_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Format data.
  data <- FormatData(
    data,
    arm_name = NULL,
    cens_after_last = cens_after_last,
    strata = strata,
    weights = weights
  )
  
  # Truncation time.
  if (is.null(tau)) {
    max_t <- NULL
    tau <- data %>% 
      dplyr::summarise(max_t = max(time)) %>% 
      dplyr::pull(max_t) %>% min()
  }

  # Calculate AUC.
  data$arm <- 0
  out <- SingleArmStratAUC(
    data = data,
    tau = tau,
    alpha = alpha,
    boot = boot,
    reps = reps
  )
  
  # Output.
  return(out)
}


# -----------------------------------------------------------------------------
# Two-sample.
# -----------------------------------------------------------------------------


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
#' @param data Data.frame. 
#' @param alpha Type I error level.
#' @param arm_name Name of column containing treatment arm. Must be coded as 1
#'   for treatment, 0 for reference.
#' @param boot Logical, construct bootstrap confidence intervals (slow)?
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
#' @param idx_name Name of column containing a unique subject index.
#' @param perm Logical, perform permutation test (slow)?
#' @param reps Number of replicates for bootstrap/permutation inference.
#' @param status_name Name of column containing the status. Must be coded as 0 for
#'   censoring, 1 for event, 2 for death. Each subject should have an 
#'   observation-terminating event, either censoring or death.
#' @param strata Optional stratification factor. Should not be provided if a
#'   covariate matrix is provided.
#' @param tau Numeric truncation time.
#' @param time_name Name of column containing the observation time.
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Object of class \code{CompStratAUCs} or \code{CompAugAUCs} with these slots:
#' \itemize{
#'   \item `@Areas`: The AUC for each arm.
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals.
#'   \item `@Curves`: Mean cumulative count curve for each arm; averaged across strata
#'     if present.
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics.
#'   \item `@Weights`: Per-stratum weights and AUCs.
#' }
#' @export
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
#'   data,
#'   tau = 2,
#'   boot = TRUE,
#'   perm = TRUE,
#'   reps = 25,
#'   alpha = 0.05
#' )
#' show(aucs)
#' }
CompareAUCs <- function(
  data,
  alpha = 0.05,
  arm_name = "arm",
  boot = FALSE,
  cens_after_last = TRUE,
  covars = NULL,
  idx_name = "idx",
  perm = FALSE,
  reps = 2000,
  status_name = "status",
  strata = NULL,
  tau = NULL,
  time_name = "time",
  weights = NULL
) {
  
  # Check that covariates and strata are not both provided.
  if (!is.null(covars) & !is.null(strata)) {
    msg <- paste0(
      "If adjustment for both strata and covariates is needed,\n",
      "include strata indicators within the covariates.")
    stop(msg)
  }
  
  # Rename columns as necessary.
  arm <- idx <- status <- time <- NULL
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      idx = {{idx_name}},
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    dplyr::select(arm, idx, status, time)
  
  # Check that no patient are in both arms.
  arm1_idx <- unique(data$idx[data$arm == 1])
  arm0_idx <- unique(data$idx[data$arm == 0])
  repeated_indices <- length(intersect(arm1_idx, arm0_idx))
  if (repeated_indices > 0) {
    msg <- paste0(
      "Patients in different arms with the same identifier (idx)\n",
      "were detected. Each patient should have a unique identifier."
    )
    stop(msg)
  }
  
  # Format data.
  data <- FormatData(
    data,
    cens_after_last = cens_after_last,
    covars = covars,
    strata = strata,
    weights = weights
  )
  
  # Truncation time.
  if (is.null(tau)) {
    max_t <- NULL
    tau <- data %>% 
      dplyr::group_by(arm) %>% 
      dplyr::summarise(max_t = max(time)) %>% 
      dplyr::pull(max_t) %>% min()
  }
  
  # Select analysis method.
  if (!is.null(covars)) {
    Analysis <- CompareAugAUCs
  } else {
    Analysis <- CompareStratAUCs
  }
  out <- Analysis(
    data = data,
    tau = tau,
    alpha = alpha,
    boot = boot,
    perm = perm,
    reps = reps
  )
  
  # Output.
  return(out)
}

