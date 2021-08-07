
#' Simulate Data for a Single Subject.
#' 
#' @param idx Subject index.
#' @param censoring_rate Rate of censoring. 
#' @param event_rate Rate of events.
#' @param death_rate Rate of terminal events. 
#' @param tau Truncation time.
#' @importFrom stats rexp
#' @return Recurrent event data for a single subject.

SimSubj <- function(
  idx,
  censoring_rate,
  event_rate,
  death_rate,
  tau
) {
  
  # Censoring and death times.
  censor <- ifelse(
    censoring_rate > 0,
    min(rexp(n = 1, rate = censoring_rate), tau),
    tau
  )
  death <- ifelse(
    death_rate > 0,
    rexp(n = 1, rate = death_rate),
    Inf
  )
  final_status <- ifelse(death <= censor, 2, 0)
  obs_time <- min(censor, death)
  
  # Simulate event times.
  events <- c()
  follow_up <- 0
  while (follow_up < obs_time) {
    
    gap <- rexp(n = 1, rate = event_rate)
    follow_up <- follow_up + gap
    
    # If follow-up at event does not exceed minimum of censoring and death,
    # then append the event.
    if (follow_up < obs_time) {
      events <- c(events, follow_up)
    }
  }
  
  # Data.
  out <- data.frame(
    "idx" = idx,
    "time" = c(events, obs_time),
    "status" = c(rep(1, length(events)), final_status)
  )
  return(out)
}


#' Simulate Recurrent Events Data
#' 
#' Simulate recurrent events data using exponential censoring, gap, and death 
#' times. Status is coded as 0 for censoring, 1 for event, 2 for death. 
#' \itemize{
#'   \item The subject-specific event rate may depend on the treatment arm,
#'   a single standard normal covariate, and a binary stratification factor.
#'   \item The subject-specific event rate is calculated as base_event_rate *
#'   exp(arm x treatment_effect + covar x covar_effect + strata x strata_effect).
#' }
#'
#' @param base_event_rate Baseline arrival rate for recurrent events.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param covar_effect Log rate ratio for the covariate, which follows a 
#'   standard normal distribution.
#' @param death_rate Arrival rate for the terminal event.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n0 Subjects in the reference arm.
#' @param n1 Subjects in the treatment arm.
#' @param strata_effect Log rate ratio for membership to the stratum.
#' @param strata_prob Probability of the binary stratification factor.
#' @param tau Truncation time.
#' @param treatment_effect Log rate ratio for treatment.
#' @return Data.frame, containing:
#' \itemize{
#'   \item Subject 'idx', treatment 'arm', a baseline covariate 'covar', 
#'     an indicator 'strata' of stratum membership.
#'  \item 
#'  \item 'time' and 'status' of the event: 0 for censoring, 1 for an 
#'     event, 2 for death.
#' }
#' 
#' @importFrom dplyr "%>%" 
#' @export

GenData <- function(
  base_event_rate = 1.0,
  censoring_rate = 0.25,
  covar_effect = log(1.25),
  death_rate = 0.25,
  frailty_variance = 0.0,
  min_event_rate = 0.05,
  n0 = 100,
  n1 = 100,
  strata_effect = log(1.25),
  strata_prob = 0.2,
  treatment_effect = log(0.5),
  tau = 4
){
  
  # Total subjects.
  n <- n1 + n0
  
  # Generate covariate data.
  covar <- data.frame(
    idx = seq_len(n),
    arm = c(rep(1, n1), rep(0, n0)),
    covar = stats::rnorm(n)
  )
  
  # Strata.
  if (strata_prob > 0) {
    covar$strata <- stats::rbinom(n = n, size = 1, prob = strata_prob)
  } else {
    covar$strata <- 0
  }
  
  # Calculate subject-specific event rate:
  covar$true_event_rate <- base_event_rate * 
    exp(covar$arm * treatment_effect + 
        covar$covar * covar_effect +
        covar$strat * strata_effect
      )
  covar$true_event_rate <- sapply(
    covar$true_event_rate, 
    function(x) {max(x, min_event_rate)}
  )
  
  # Frailty.
  if (frailty_variance > 0) {
    theta <- 1 / frailty_variance
    covar$frailty <- stats::rgamma(n = n, shape = theta, rate = theta)
  } else {
    covar$frailty <- 1
  }
  covar$true_event_rate <- covar$frailty * covar$true_event_rate 
  covar$true_death_rate <- covar$frailty * death_rate
  
  # Generate event-date.
  idx <- NULL
  true_death_rate <- NULL
  true_event_rate <- NULL
  data <- covar %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      SimSubj(
        idx = idx, 
        censoring_rate = censoring_rate, 
        event_rate = true_event_rate, 
        death_rate = true_death_rate, 
        tau = tau
      ),
      .groups = "drop"
    ) %>% as.data.frame()
  
  # Merge in covariates.
  data <- merge(
    x = data,
    y = covar,
    by = "idx"
  )
  return(data)
}