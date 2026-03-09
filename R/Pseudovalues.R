# Purpose: Calculations of the AUC.
# Updated: 2026-03-08


#' Generate Pseudovalues
#'
#' @param data Data.frame.
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event?
#' @param idx_name Name of column containing a unique subject index.
#' @param status_name Name of column containing the status. Must be coded as 0 for
#'   censoring, 1 for event, 2 for death. Each subject should have an
#'   observation-terminating event, either censoring or death.
#' @param tau Numeric truncation time.
#' @param time_name Name of column containing the observation time.
#' @param type Character. \code{"MCF"} for pseudo-values of the MCF at \code{tau},
#'   or \code{"AUC"} for pseudo-values of the area under the MCF up to \code{tau}.
#'   Both give one row per subject.
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Data.frame with \code{idx}, \code{psi}, and \code{pseudo} (one row per subject).
#' @export
GenPseudo <- function(
    data,
    cens_after_last = TRUE,
    idx_name = "idx",
    status_name = "status",
    tau = NULL,
    time_name = "time",
    type = c("MCF", "AUC"),
    weights = NULL
) {

  type <- match.arg(type)

  # Rename columns as necessary.
  idx <- orig_idx <- time <- status <- NULL
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
    weights = weights,
    keep_orig_idx = TRUE
  )
  idx_map <- data %>%
    dplyr::select(idx, orig_idx) %>%
    unique()

  # Truncation time.
  if (is.null(tau)) {
    max_t <- NULL
    tau <- data %>%
      dplyr::summarise(max_t = max(time)) %>%
      dplyr::pull(max_t) %>%
      min()
  }

  # Tabulate MCF (needed for both types).
  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    weights = data$weights,
    calc_var = FALSE
  )

  # Calculate influence function.
  out <- InfluenceFunction(
    data = data,
    tau = tau,
    type = type,
    mcf = mcf,
    weights = data$weights
  )

  # Calculate parameter.
  if (type == "AUC") {
    auc <- SingleArmAUC(data = data, tau = tau)
    param <- auc@MargAreas$area[1]
  } else {
    tau_eff <- max(mcf$time[mcf$time <= tau], na.rm = TRUE)
    if (!is.finite(tau_eff)) {
      tau_eff <- min(mcf$time)
    }
    param <- mcf$mcf[mcf$time == tau_eff][1]
  }

  # Construct pseudovalues.
  out$pseudo <- param + out$psi

  # Format output.
  out <- out %>%
    dplyr::inner_join(idx_map, by = "idx") %>%
    dplyr::select(-idx) %>%
    dplyr::relocate(orig_idx) %>%
    dplyr::rename(idx = orig_idx)
  return(out)
}
