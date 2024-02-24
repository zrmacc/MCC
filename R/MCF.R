# Purpose: R wrapper for CalcMCF function.
# Updated: 2024-02-21

#' Calculate Mean Cumulative Function
#' 
#' Tabulates the mean cumulative function. See equation 2.1 of 
#'  <doi:10.1111/j.0006-341X.2000.00554.x>.
#'  
#' @param idx Unique subject index. 
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#' @param time Observation time.
#' @param calc_var Calculate variance of the MCF?
#' @param weights Jump weights.
#' @return Data.frame with these columns:
#' \itemize{
#'    \item `times`, distinct observation times.
#'    \item `censor`, number of censorings.
#'    \item `death`, number of deaths.
#'    \item `event`, number of events.
#'    \item `haz`, instantaneous hazard (of death).
#'    \item `surv`, survival probability.
#'    \item `event_rate`, instantaneous event rate.
#'    \item `mcf`, mean cumulative function. 
#'    \item `se_mcf`, standard error of the MCF.
#' }
#' @export
CalcMCF <- function(
  idx,
  status,
  time,
  calc_var = TRUE,
  weights = NULL
){
  
  # Weights.
  n <- length(idx)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  
  # Call MCF.
  out <- CalcMCFCpp(
    idx = idx,
    status = status,
    time = time,
    weights = weights,
    calc_var = calc_var
  )
  return(out)
}


