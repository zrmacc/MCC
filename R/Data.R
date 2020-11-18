
#' Simulate Data for a Single Subject.
#' 
#' @param idx Subject index.
#' @param censoring_rate Rate of censoring. 
#' @param event_rate Rate of events.
#' @param death_rate Rate of terminal events. 
#' @param tau Truncation time.

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
#' @param n Numbers of subjects.
#' @param censoring_rate Rate of censoring. 
#' @param event_rate Rate of events.
#' @param death_rate Rate of terminal events. 
#' @param tau Truncation time.
#' @importFrom stats rexp
#' @export
#' @return Data.frame containing:
#' \itemize{
#'   \item Subject 'idx'.
#'   \item Observation 'time'.
#'   \item Observation 'status', 0 for censoring, 1 for event, 2 for death.
#' }

GenData <- function(
  n, 
  censoring_rate, 
  event_rate, 
  death_rate, 
  tau
) {
  
  # Wrap per-subject data generation.
  aux <- function(i) {
    SimSubj(
      idx = i,
      censoring_rate = censoring_rate,
      event_rate = event_rate,
      death_rate = death_rate,
      tau = tau
    )
  }
  
  # Generate data.
  data <- lapply(
    seq_len(n),
    aux
  )
  data <- do.call(rbind, data)
  
  # Return.
  return(data)
}


# -----------------------------------------------------------------------------

#' Sample Recurrent Event Data
#'
#' Sample recurrent events data. Events occurred in the treatment arm according
#' to a Poisson process with rate 1.00, and in the reference arm according to a
#' Poisson process with rate 1.35. Data were truncated at time 10. Death and
#' censoring occurred independently according to an exponential distribution
#' with rate 0.25 in each arm. The number of subjects in each arm is 100.
#'
#' @docType data
#' @usage data(mcc_data)
#' @format A data.frame containing three fields:
#' \describe{
#'   \item{idx}{Unique subject identifier.}
#'   \item{time}{Observation time between 0 and 10.}
#'   \item{arm}{Treatment arm, 0 for reference, 1 for treatment.}
#'   \item{status}{Status indicator, 0 for censoring, 1 for an event, 2 for death.}
#' }
"mcc_data"