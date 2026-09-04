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
#' @param jump_weights Optional column of jump weights, controlling the size of the jump
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
    jump_weights = NULL
) {

  prepared <- .PreparePseudo(
    data = data,
    cens_after_last = cens_after_last,
    idx_name = idx_name,
    status_name = status_name,
    tau = tau,
    time_name = time_name,
    type = type,
    jump_weights = jump_weights
  )
  return(prepared$pseudo)
}


#' Prepare pseudo-values and their process inputs
#'
#' @return Internal list used by \code{GenPseudo()} and \code{PseudoReg()}.
#' @noRd
.PreparePseudo <- function(
    data,
    cens_after_last = TRUE,
    idx_name = "idx",
    status_name = "status",
    tau = NULL,
    time_name = "time",
    type = c("MCF", "AUC"),
    jump_weights = NULL
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
    jump_weights = jump_weights,
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
    jump_weights = data$jump_weights,
    calc_var = FALSE
  )

  # Calculate influence function.
  out <- InfluenceFunction(
    data = data,
    tau = tau,
    type = type,
    mcf = mcf,
    jump_weights = data$jump_weights
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
  return(list(
    pseudo = out,
    process_data = data,
    idx_map = idx_map,
    mcf = mcf,
    parameter = param,
    tau = tau,
    tau_effective = if (type == "MCF") tau_eff else tau,
    type = type
  ))
}
