# Purpose: Functions to perform inference on the difference or ratio 
# in AUCs between two mean cumulative count curves. 

# -----------------------------------------------------------------------------
# Area under the mean cumulative function.
# -----------------------------------------------------------------------------

#' Find Area Under the Mean Cumulative Count Curve. 
#'
#' @param times MCF times.
#' @param values MCF values
#' @param tau Truncation time.
#' @importFrom stats integrate stepfun 
#' @return Numeric area under the curve. 

FindAUC <- function(times, values, tau) {
  g <- stepfun(
    x = times,
    y = c(0, values)
  )
  area <- integrate(f = g, lower = 0, upper = tau, subdivisions = 2e3)
  return(area$value)
}

# -----------------------------------------------------------------------------
# Variance calculation for area under the mean cumulative function.
# -----------------------------------------------------------------------------

#' Calculate Influence Function Contribution of a Single Subject.
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times for a single subject.
#' @param status Status indicator for a single subject.
#' @param tau Truncation time.
#' @return Numeric influence function for a single observation.

CalcPsi.AUC.i <- function(
  mcf,
  time,
  status,
  tau
) { 

  # Subject-specific martingales.
  out <- CalcMartigales.i(
    mcf = mcf,
    time = time,
    status = status
  )
  
  # Truncate.
  out <- out[out$time <= tau, ]
  mcf <- mcf[mcf$time <= tau, ]
  
  # Calculate variance.
  int_1 <- sum((tau - out$time) * mcf$surv / mcf$prop_at_risk * out$dM_event)
  int_2 <- sum((tau - out$time) * mcf$mcf  / mcf$prop_at_risk * out$dM_death)
  psi <- int_1 + int_2
  
  # Output.
  return(psi)
}


#' Calculate Variance of the Area Under the Mean Cumulative Function
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times.
#' @param status Status indicators.
#' @param idx Unique subject index. 
#' @param tau Truncation time.
#' @importFrom dplyr "%>%" group_by select summarise
#' @return Numeric influence function contributions of each subject.

CalcPsi.AUC <- function(
  mcf,
  time,
  status,
  idx,
  tau
) {
  
  subj <- data.frame(
    "time" = time,
    "status" = status,
    "idx" = idx
  )
  
  # Generate data.frame with each subject's influence function.
  out <- subj %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      "psi" = CalcPsi.AUC.i(mcf = mcf, time = time, status = status, tau = tau),
      .groups = "drop" ) %>%
    as.data.frame

  # Output.
  return(out)
}

# -----------------------------------------------------------------------------

#' Calculate AUC
#' 
#' Used to find the AU MCF for each stratum.
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period. 
#' @param idx Unique subject index. 
#' @param tau Truncation time. 
#' @param calc_var Calculate analytical variance? 
#' @return Data.frame containing:
#' \itemize{
#'   \item Truncation time 'tau' and estimated 'area'.
#'   \item Variance 'var_area' and standard error 'se_area', if requested.
#' }

AUC.Area <- function(
  time, 
  status,
  idx,
  tau,
  calc_var = TRUE
) { 
  
  # Fit MCF.
  fit <- CalcMCF(
    time = time,
    status = status,
    idx = idx,
    calc_var = calc_var
  )
  
  # Find Area.
  area <- FindAUC(
    times = fit$time,
    values = fit$mcf,
    tau = tau
  )
  
  # Output.
  n <- length(unique(idx))
  out <- data.frame(
    "n" = n,
    "tau" = tau,
    "area" = area
  )
  
  if (calc_var) {
    
    # Calculate influence contributions.
    psi <- CalcPsi.AUC(
      mcf = fit,
      time = time,
      status = status,
      idx = idx,
      tau = tau
    )
    
    # Find variance of area.
    out$var_area <- mean(psi$psi^2)
    out$se_area <- sqrt(out$var_area / n)
  }
  
  return(out)
}