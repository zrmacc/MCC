# Purpose: Calculations of the AUC.
# Updated: 2020-12-17

# -----------------------------------------------------------------------------
# Area under the mean cumulative function.
# -----------------------------------------------------------------------------

#' Find Area Under the Curve
#'
#' @param times MCF times.
#' @param values MCF values
#' @param tau Truncation time.
#' @importFrom stats integrate stepfun 
#' @return Numeric area under the curve. 

AUC <- function(times, values, tau) {
  g <- stepfun(
    x = times,
    y = c(0, values)
  )
  area <- integrate(f = g, lower = 0, upper = tau, subdivisions = 1e4)
  return(area$value)
}

# -----------------------------------------------------------------------------
# Influence function for AUC.
# -----------------------------------------------------------------------------

#' Calculate Influence Function Contribution of a Single Subject.
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times for a single subject.
#' @param status Status indicator for a single subject.
#' @param tau Truncation time.
#' @return Numeric influence function for a single observation.

PsiAUCi <- function(
  mcf,
  time,
  status,
  tau
) { 

  # Subject-specific martingales.
  out <- Martingales(
    mcf = mcf,
    time = time,
    status = status
  )
  
  # Truncate.
  out <- out[out$time <= tau, ]
  mcf <- mcf[mcf$time <= tau, ]
  
  # First integral.
  int_1 <- sum((tau - out$time) * mcf$surv / mcf$prop_at_risk * out$dM_event)
  
  # Second integral.
  out$nu <- sum((tau - out$time) * mcf$surv * mcf$event_rate) - 
    cumsum((tau - out$time) * mcf$surv * mcf$event_rate)
  int_2 <- sum(out$nu / mcf$prop_at_risk * out$dM_death)
  
  # Influence function.
  psi <- int_1 - int_2

  # Output.
  return(psi)
}

