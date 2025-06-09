# Purpose: Calculations of the AUC.
# Updated: 2025-06-08

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
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Data.frame with the influence and pseudo-value for each unique observation.
#' @export
GenPseudo <- function(
    data,
    cens_after_last = TRUE,
    idx_name = "idx",
    status_name = "status",
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
  auc <- SingleArmAUC(data = data, tau = tau)
  param <- auc@MargAreas$area[1]
  
  # Tabulate MCF.
  mcf <- CalcMCF(
    idx = data$idx, 
    status = data$status, 
    time = data$time, 
    weights = data$weights, 
    calc_var = FALSE
  )
  
  # Calculate influence.
  data_tau <- data %>% dplyr::filter(time <= tau)
  mcf_tau <- mcf %>% dplyr::filter(time <= tau)
  
  out <- PsiAUC(
    event_rate = mcf_tau$weighted_event_rate,
    idx = data_tau$idx,
    haz = mcf_tau$haz,
    nar = mcf_tau$nar,
    status = data_tau$status,
    surv = mcf_tau$surv,
    tau = tau,
    time = data_tau$time,
    weights = data_tau$weights
  )
  
  # Output.
  out$pseudo <- param + out$psi
  return(out)
}
