# Purpose: Estimate the mean cumulative function using the method of Ghosh and Lin.

# -----------------------------------------------------------------------------

#' Calculate Mean Cumulative Function
#' 
#' Tabulates the mean cumulative function. See equation 2.1 of 
#' <doi:10.1111/j.0006-341X.2000.00554.x>.
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#' @param idx Unique subject index. 
#' @export 
#' @return Data.frame with these columns:
#' \itemize{
#'   \item `times`, distinct observation times.
#'   \item `censor`, number of censorings.
#'   \item `death`, number of deaths.
#'   \item `event`, number of events.
#'   \item `haz`, instantaneous hazard (of death).
#'   \item `surv`, survival probability.
#'   \item `event_rate`, instantaneous event rate.
#'   \item `mcf`, mean cumulative function. 
#' }

CalcMCF <- function(
  time,
  status,
  idx
) {
  
  # Format data.
  mcf <- TabulateEvents(time, status, idx)
  
  # Kaplan-Meier estimator.
  mcf$haz <- mcf$death / mcf$nar
  mcf$surv <- cumprod(1 - mcf$haz)
  
  # Nelson-Aalen estimator of event rate.
  mcf$event_rate <- mcf$event / mcf$nar
  
  # Estimate marginal cumulative function.
  mcf$mcf <- cumsum(mcf$surv * mcf$event_rate)
  
  # Calculate standard error.
  n <- length(unique(idx))
  mcf$var_mcf <- CalcVarMCF(
    mcf = mcf,
    time = time,
    status = status,
    idx = idx
  )
  mcf$se_mcf <- sqrt(mcf$var_mcf / n)
  
  # Output.
  return(mcf)
}


# -----------------------------------------------------------------------------

#' Tabulate Events
#' 
#' Tabulate the number at risk and the number of events at each unique
#' observation time.
#' 
#' @param time Event time.
#' @param status Status, coded as 1 for the recurrent event, 0 otherwise.
#' @param idx Subject index. 
#' @return Data.frame with these columns:
#' \itemize{
#'   \item `time` distinct observation times.
#'   \item `censor` censorings.
#'   \item `deaths` deaths.
#'   \item `events` events.
#' }

TabulateEvents <- function(time, status, idx){
  
  # Censoring times.
  censor_times <- sort(time[status == 0])
  
  # Event times.
  event_times <- sort(time[status == 1])
  
  # Death times.
  death_times <- sort(time[status == 2])
  
  # Unique observation times.
  out <- data.frame("time" = sort(unique(time)))
  
  # Number initially at risk. 
  init_nar <- length(unique(idx))
  
  # Number of censorings.
  out$censor <- 0
  if (length(censor_times) > 0) {
    out$censor[out$time %in% censor_times] <- 1
  }
  
  # Number of deaths.
  out$death <- 0
  if (length(death_times) > 0) {
    out$death[out$time %in% death_times] <- 1
  }
  
  # Number events.
  out$events <- 0
  if (length(event_times) > 0) {
    out$events[out$time %in% event_times] <- 1
  }
  
  # Number at risk.
  out$nar <- sapply(out$time, FUN = function(time) {
    return(init_nar - sum(death_times < time) - sum(censor_times < time))
  })
  out$prop_at_risk <- out$nar / init_nar
  
  # Output.
  return(out)
}


# -----------------------------------------------------------------------------

#' Calculate Martingales for a Single Subject
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times for a single subject.
#' @param status Status indicator for a single subject.
#' @return Data.frame containing:
#' \itemize{
#'   \item Distinct observation times from `mcf`.
#'   \item 'censor', 'death', and 'event' indicators.
#'   \item 'nar', an at-risk indicator.
#'   \item 'dM_event' and 'dM_death', increments of the event and 
#'     death martingales.
#' }

CalcMartigales.i <- function(
  mcf,
  time,
  status
) {
  
  # Censoring time.
  censor_times <- time[status == 0]
  
  # Event times.
  event_times <- sort(time[status == 1])
  
  # Death times.
  death_times <- sort(time[status == 2])
  
  # Subject specific events.
  out <- data.frame("time" = mcf$time)
  
  # Censoring indicator.
  out$censor <- 0
  if (length(censor_times) > 0) {
    out$censor[out$time == censor_times] <- 1
  }
  
  # Death indicator.
  out$death <- 0
  if (length(death_times) > 0) {
    out$death[out$time == death_times] <- 1
  }
  
  # Event indicator.
  out$event <- 0
  if (length(event_times) > 0) {
    out$event[out$time %in% event_times] <- 1
  }
  
  # Subject specific at risk.
  out$nar <- sapply(out$time, FUN = function(time){
    return(1 - sum(death_times < time) - sum(censor_times < time))
  })
  
  # Calculate martingale differences
  out$dM_event <- out$event - out$nar * mcf$event_rate
  out$dM_death <- out$death  - out$nar * mcf$haz
  
  # Output.
  return(out)
}

# -----------------------------------------------------------------------------

#' Calculate MCF Variance Contribution for a Single Subject
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times for a single subject.
#' @param status Status indicator for a single subject.
#' @return Numeric vector containing the squared influence function for a single
#'   subject at each distinct event event.

CalcVarMCF.i <- function(
  mcf,
  time,
  status
) { 
  
  # Subject-specific martingales.
  out <- CalcMartigales.i(
    mcf = mcf,
    time = time,
    status = status
  )
  
  # Calculate squared influence function.
  # Integrals correspond to those in Ghosh and Lin (2.2).
  out$int_1 <- cumsum(mcf$surv / mcf$prop_at_risk * out$dM_event)
  out$int_2 <- mcf$mcf * cumsum(out$dM_death / mcf$prop_at_risk)
  out$int_3 <- cumsum(mcf$mcf * out$dM_death / mcf$prop_at_risk)
  out$psi <- out$int_1 - out$int_2 + out$int_3
  out$psi2 <- (out$psi)^2
  return(out$psi2)
}


#' Calculate variance of MCF.
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times.
#' @param status Status indicators.
#' @param idx Unique subject index. 
#' @return Numeric vector of estimated variances at each distinct observation
#'   time in `mcf`.

CalcVarMCF <- function(
  mcf,
  time,
  status,
  idx
) {
  
  # Partition data by subject.
  subj <- data.frame(
    "time" = time,
    "status" = status,
    "idx" = idx
  )
  subj_split <- split(subj, subj$idx)
  
  # Calculate psi squared; Ghosh and Lin (2.2).
  psi2 <- lapply(subj_split, FUN = function(df) {
    return(CalcVarMCF.i(mcf = mcf, time = df$time, status = df$status))
  })
  var <- Reduce("+", psi2) / length(psi2)
  
  # Output.
  return(var)
}

# -----------------------------------------------------------------------------

#' Calculate Average MCF Curve
#' 
#' Calculates the weighted average of MCF curves.
#' 
#' @param curve_list List of tabulated MCFs as returned by \code{\link{CalcMCF}}.
#' @param weights Numeric vector of weights.
#' @return Data.frame containing `Time` and the average MCF `Avg_MCF`.

AvgMCF <- function (curve_list, weights) {
  
  # Extract event times.
  time <- lapply(curve_list, function(x){x$time})
  time <- do.call(c, time)
  time <- sort(unique(time))
  
  # Extract mcfs evaluated on time.
  aux <- function(x){
    g <- stepfun(x$time, c(0, x$mcf), right = FALSE)
    return(g(time))
  }
  mcfs <- lapply(curve_list, aux)
  mcfs <- do.call(cbind, mcfs)
  avg_mcf <- mcfs %*% weights
  
  # Standard errors.
  aux <- function(x){
    g <- stepfun(x$time, c(0, x$var_mcf), right = FALSE)
    return(g(time))
  }
  vars <- lapply(curve_list, aux)
  vars <- do.call(cbind, vars)
  avg_var <- vars %*% (weights^2)
  
  # Output table.
  out <- data.frame(
    'time' = time,
    'mcf' = avg_mcf,
    'var_mcf' = avg_var
  )
  out$se_mcf <- sqrt(out$var)
  return(out)
}

