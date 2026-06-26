# Purpose: Main function for multivariate AUC comparison.
# Updated: 2026-06-16

#' Inference on Multivariate Areas Under the Cumulative Count Curves
#'
#' Confidence intervals and p-values for the multivariate difference in areas
#' under the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0). By default returns the K-dimensional contrast vector
#' \eqn{\hat{\Delta}(\tau)}; optional \code{process_weights} provide a
#' secondary scalar combination \eqn{w^\top\hat{\Delta}(\tau)} in the
#' \code{Weighted} slot.
#'
#' @section Terminal events:
#' Data are in long format with an \code{event_type} column on every row.
#' Recurrent events (\code{status == 1}) must have a non-missing event type.
#' For each event type \eqn{k} for which a subject is at risk (see section
#' \strong{Risk sets}), follow-up must end with a single observation-terminating
#' event (\code{status} is 0 or 2). A terminal row with missing
#' \code{event_type} ends follow-up for all event types; a terminal row with
#' \code{event_type = k} ends follow-up for process \eqn{k} only. If a patient
#' has a terminal event for some event types but not others, then by default
#' (\code{cens_after_last = TRUE}) a censoring event is added immediately after
#' the last recurrent event for each remaining at-risk event type. Set
#' \code{cens_after_last = FALSE} to leave such patients at risk indefinitely
#' (with a warning).
#'
#' @section Risk sets:
#' A subject is at risk for process \eqn{k} if they have at least one row with
#' non-missing \code{event_type = k} (recurrent or typed terminal). Subjects
#' with no type-\eqn{k} rows are excluded from the analysis of \eqn{k}. A global
#' terminal row does not place a subject in the risk set for processes without
#' typed rows. To include a subject at risk for \eqn{k} with zero recurrent
#' type-\eqn{k} events, add a typed censoring or death row with
#' \code{event_type = k}. Drop event rows that should not count before calling
#' this function.
#'
#' @section Inference:
#' When \code{reps = NULL} (default), Wald confidence intervals are based on
#' the influence-function covariance. These are marginal (per event type), not
#' simultaneous over \eqn{\hat{\Delta}}. When \code{reps} is a positive
#' integer, subject-level bootstrap intervals and p-values are appended.
#'
#' @param data Data.frame in long format with columns for arm, subject index,
#'   time, status, and event type.
#' @param alpha Type I error level.
#' @param arm_name Name of treatment arm column (1 = treatment, 0 = reference).
#' @param cens_after_last Censor subjects without a terminal row after their
#'   last recurrent event?
#' @param covars Optional subject-level covariates for augmentation adjustment.
#' @param event_type_name Name of event type column.
#' @param idx_name Name of subject index column.
#' @param reps Number of bootstrap replicates; \code{NULL} for asymptotic only.
#' @param status_name Name of status column (0 = censor, 1 = event, 2 = death).
#' @param tau Truncation time; default is the minimum of the per-arm maximum
#'   observation times.
#' @param time_name Name of time column.
#' @param process_weights Optional numeric vector of length \eqn{K} defining
#'   contrast weights for the secondary scalar estimand
#'   \eqn{w^\top\hat{\Delta}(\tau)} in the \code{Weighted} slot.
#' @return Object of class \code{CompMvAugAUCs}.
#' @seealso \code{\link{CompareAUCs}} for univariate analysis.
#' @export
#' @examples
#' \donttest{
#' # K = 1 (backward compatible with univariate data).
#' covariates <- data.frame(arm = c(rep(1, 50), rep(0, 50)))
#' data <- GenData(beta_event = log(0.5), covariates = covariates)
#' mv <- CompareMvAUCs(data, tau = 2)
#' show(mv)
#' }
CompareMvAUCs <- function(
    data,
    alpha = 0.05,
    arm_name = "arm",
    cens_after_last = TRUE,
    covars = NULL,
    event_type_name = "event_type",
    idx_name = "idx",
    reps = NULL,
    status_name = "status",
    tau = NULL,
    time_name = "time",
    process_weights = NULL
) {

  data <- data %>%
    dplyr::rename(
      arm = dplyr::all_of(arm_name),
      idx = dplyr::all_of(idx_name),
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
    )

  arm1_idx <- unique(data[["idx"]][data[["arm"]] == 1])
  arm0_idx <- unique(data[["idx"]][data[["arm"]] == 0])
  if (length(intersect(arm1_idx, arm0_idx)) > 0) {
    stop(
      "Patients in different arms with the same identifier (idx) ",
      "were detected. Each patient should have a unique identifier."
    )
  }

  fmt <- FormatMvData(
    data = data,
    arm_name = "arm",
    cens_after_last = cens_after_last,
    covars = covars,
    event_type_name = event_type_name,
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
  )

  data <- fmt$data
  event_types <- fmt$event_types
  cens_after_last <- if (fmt$univariate_path) {
    cens_after_last
  } else {
    fmt$cens_after_last
  }

  covar_cols <- NULL
  if (!is.null(covars)) {
    covar_cols <- setdiff(
      names(data),
      c(
        "arm", "idx", "status", "time", "jump_weights", "strata",
        "event_type", "orig_idx"
      )
    )
  }

  if (is.null(tau)) {
    tau <- data %>%
      dplyr::group_by(arm) %>%
      dplyr::summarise(max_t = max(time), .groups = "drop") %>%
      dplyr::pull(max_t) %>%
      min()
  }

  obs <- CalcMvAugAUC(
    data = data,
    event_types = event_types,
    tau = tau,
    alpha = alpha,
    cens_after_last = cens_after_last,
    covar_cols = covar_cols,
    process_weights = process_weights,
    return_areas = TRUE
  )

  obs_stats <- obs$contrasts
  res <- .AddMvResamplingResults(
    obs_stats = obs_stats,
    data = data,
    event_types = event_types,
    tau = tau,
    alpha = alpha,
    reps = reps,
    cens_after_last = cens_after_last,
    covar_cols = covar_cols,
    process_weights = process_weights
  )

  weighted <- obs$weighted
  if (is.null(process_weights)) {
    weighted <- data.frame()
  }

  methods::new(
    Class = "CompMvAugAUCs",
    Areas = obs$marg_areas,
    CIs = res$cis,
    CovDelta = obs$cov_delta,
    CovTheta = obs$cov_theta,
    MCF = obs$mcf,
    Pvals = res$pvals,
    Reps = res$sim_reps,
    Weighted = weighted
  )
}
