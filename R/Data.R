
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
  censor <- rexp(n = 1, rate = censoring_rate)
  death <- rexp(n = 1, rate = death_rate)
  min_censor_death <- min(censor, death)
  final_status <- ifelse(min_censor_death == death, 2, 0)
  
  # Add truncation.
  min_censor_death_tau <- min(min_censor_death, tau)
  final_status <- ifelse(min_censor_death_tau == tau, 0, final_status)
  
  # Simulate event times.
  events <- c()
  follow_up <- 0
  while (follow_up < min_censor_death) {
    
    gap <- rexp(n = 1, rate = event_rate)
    follow_up <- follow_up + gap
    
    # If follow-up at event does not exceed minimum of censoring and death,
    # then append the event.
    if (follow_up < min_censor_death_tau) {
      events <- c(events, follow_up)
    }
  }
  
  # Data.
  out <- data.frame(
    "idx" = idx,
    "time" = c(events, min_censor_death_tau),
    "status" = c(rep(1, length(events)), final_status)
  )
  return(out)
}


#' Simulate Recurrent Events Data
#' 
#' Simulate recurrent events data using exponential censoring, gap, and death 
#' times. Status is coded as 0 for censoring, 1 for event, 2 for death. 
#'
#' @param n1 Subjects in the treatment arm.
#' @param n0 Subjects in the reference arm.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param base_event_rate Baseline arrival rate for recurrent events.
#' @param death_rate Arrival rate for the terminal event.
#' @param treatment_effect Log rate ratio for treatment.
#' @param pi Probability of membership to the high-risk stratum.
#' @param risk_effect Log rate ratio for membership to the high risk stratum.
#' @param covar_effect Log rate ratio for the covariate, which follows a 
#'   standard normal distribution.
#' @param min_event_rate Minimum subject-specific event rate. Most be positive.
#' @param tau Truncation time.
#' @param seed Data generation seed.
#' @import dplyr
#' @importFrom stats rnorm
#' @export
#' @return Data.frame, containing:
#' \itemize{
#'   \item Subject 'idx', treatment 'arm', a baseline covariate 'covar', 
#'     an indicator 'strat' of membership to the high-risk stratum, and the
#'     subject's true underlying event rate.
#'   \item 'time' and 'status' of the event: 0 for censoring, 1 for an 
#'     event, 2 for death.
#' }

GenData <- function(
  n1 = 100,
  n0 = 100,
  censoring_rate = 0.25,
  base_event_rate = 1.0,
  death_rate = 0.25,
  treatment_effect = log(0.5),
  pi = 0.2,
  risk_effect = log(1.25),
  covar_effect = log(1.25),
  min_event_rate = 0.05,
  tau = 4,
  seed = 101
){
  
  # Set seed.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Total subjects.
  n <- n1 + n0
  
  # Generate covariate data.
  covar <- data.frame(
    "idx" = seq_len(n),
    "arm" = c(rep(1, n1), rep(0, n0)),
    "covar" = rnorm(n),
    "strat" = rbinom(n = n, size = 1, prob = pi)
  )
  
  # Calculate subject-specific event rate:
  covar$true_event_rate <- base_event_rate * 
    exp(
      covar$arm * treatment_effect + 
        covar$strat * risk_effect +
        covar$covar * covar_effect 
      )
  covar$true_event_rate <- sapply(
    covar$true_event_rate, 
    function(x) {max(x, min_event_rate)}
  )
  
  # Generate event-date.
  idx <- true_event_rate <- NULL
  data <- covar %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      SimSubj(
        idx = idx, 
        censoring_rate = censoring_rate, 
        event_rate = true_event_rate, 
        death_rate = death_rate, 
        tau = tau
      ),
      .groups = "drop"
    ) %>% as.data.frame
  
  # Merge in covariates.
  data <- merge(
    x = data,
    y = covar,
    by = "idx"
  )
  return(data)
}