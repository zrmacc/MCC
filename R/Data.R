# Purpose: Data generation.
# Updated: 2022-05-15

#' Simulate Recurrent Events Data
#' 
#' Simulate recurrent events data using exponential censoring, gap, and death 
#' times. Status is coded as 0 for censoring, 1 for event, 2 for death.
#' \itemize{
#'   \item The subject-specific death rate is calculated as frailty x base_death_rate x 
#'   exp(beta_death x covariates).
#'   \item The subject-specific event rate is calculated as frailty x base_event_rate x
#'   exp(beta_event x covariates).
#' }
#'
#' @param base_death_rate Baseline arrival rate for the terminal event.
#' @param base_event_rate Baseline arrival rate for recurrent events.
#' @param beta_death Numeric vector of log rate ratios for the death rate.
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param covariates Numeric design matrix.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_death_rate Minimum subject-specific event rate. Must be positive.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n Number of subjects. Overwritten by `nrow(covariates)` if covariates are provided.
#' @param tau Truncation time.
#' @return Data.frame, containing:
#' \itemize{
#'  \item The subject identifier `idx` (index).
#'  \item `time` and `status` of the event: 0 for censoring, 1 for an 
#'     event, 2 for death.
#'  \item The `true_death_rate`, `true_event_rate`, and `frailty`, which are subject-specific.
#' }
#' 
#' @importFrom dplyr "%>%" 
#' @export

GenData <- function(
  base_death_rate = 0.25,
  base_event_rate = 1.0,
  beta_death = NULL,
  beta_event = NULL,
  censoring_rate = 0.25,
  covariates = NULL,
  frailty_variance = 0.0,
  min_death_rate = 0.05,
  min_event_rate = 0.05,
  n = 100,
  tau = 4
){
  
  # Generate subject-specific data frame.
  if (is.null(covariates)) {
    covariates <- data.matrix(rep(1, n))
  } else {
    covariates <- data.matrix(covariates)
    n <- nrow(covariates)
  }
  df <- data.frame(idx = seq_len(n), covariates)
  
  # Calculate subject-specific event rate.
  if (is.null(beta_death)) {beta_death <- rep(0, ncol(covariates))}
  if (is.null(beta_event)) {beta_event <- rep(0, ncol(covariates))}
  df$cens_rate <- censoring_rate
  df$death_rate <- base_death_rate * exp(covariates %*% beta_death)
  df$event_rate <- base_event_rate * exp(covariates %*% beta_event)
  
  # Apply floor to death and event rates.
  df$death_rate <- pmax(df$death_rate, min_death_rate)
  df$event_rate <- pmax(df$event_rate, min_event_rate)
  
  # Apply frailty.
  if (frailty_variance > 0) {
    theta <- 1 / frailty_variance
    df$frailty <- stats::rgamma(n = n, shape = theta, rate = theta)
  } else {
    df$frailty <- 1
  }
  df$death_rate <- df$frailty * df$death_rate
  df$event_rate <- df$frailty * df$event_rate 
  
  # Generate event-date.
  data <- SimDataCpp(
    df$cens_rate, df$death_rate, df$idx, df$event_rate, tau)
  
  # Merge in covariates.
  out <- merge(
    x = data,
    y = df,
    by = "idx"
  )
  return(out)
}
