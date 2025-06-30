# Purpose: Calculations of the AUC.
# Updated: 2024-02-19

# -----------------------------------------------------------------------------
# Area under the mean cumulative function.
# -----------------------------------------------------------------------------

#' Find Area Under the Curve
#'
#' @param times Time points.
#' @param values Values.
#' @param tau Truncation time.
#' @return Numeric area under the curve.
#' @export
AUC <- function(times, values, tau) {
  
  # Extend times to include 0 and tau.
  times_ext <- c(0, times, tau)
  values_ext <- c(0, values)
  deltas <- diff(values_ext)
  if (any(deltas < 0)) {
    stop("The values supplied to AUC should be monotone increasing.")
  }
  
  # Step intervals.
  lefts <- head(times_ext, -1)
  rights <- tail(times_ext, -1)
  heights <- values_ext
  
  # Clip intervals to [0, tau].
  lefts <- pmax(lefts, 0)
  rights <- pmin(rights, tau)
  
  # Compute widths and contribution
  widths <- rights - lefts
  area <- sum(widths[widths > 0] * heights[widths > 0])
  return(area)
}


#' Calculate Variance of AUC
#' 
#' @param data Data.frame containing {idx, status, time}.
#' @param tau Truncation time.
#' @param mcf Tabulated MCF, if already computed.
#' @param return_psi Return influence function contributions?
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Numeric variance.
#' @export
VarAUC <- function(data, tau, mcf = NULL, return_psi = FALSE, weights = NULL) {
  
  if (is.null(weights)) {weights <- 1}
  data$weights <- weights
  
  if (is.null(mcf)) {
    mcf <- CalcMCF(
      idx = data$idx, 
      status = data$status, 
      time = data$time, 
      weights = data$weights, 
      calc_var = FALSE
    )
  }
  
  # Truncate.
  time <- NULL
  data_tau <- data %>% dplyr::filter(time <= tau)
  mcf_tau <- mcf %>% dplyr::filter(time <= tau)
  
  # Variance calculation.
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
  
  if (return_psi) {
    return(out)
  } else {
    return(mean(out$psi^2))
  }
}

