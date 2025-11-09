# Purpose: Calculations of the AUC.
# Updated: 2025-11-08

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
  lefts <- utils::head(times_ext, -1)
  rights <- utils::tail(times_ext, -1)
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
#' @param data Data.frame containing (idx, status, time).
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
  
  # Variance calculation.
  out <- PsiAUC(
    event_rate = mcf$weighted_event_rate,
    grid_time = mcf$time,
    idx = data$idx,
    haz = mcf$haz,
    nar = mcf$nar,
    status = data$status,
    surv = mcf$surv,
    tau = tau,
    time = data$time,
    weights = data$weights
  )
  
  if (return_psi) {
    return(out)
  } else {
    return(mean(out$psi^2))
  }
}

