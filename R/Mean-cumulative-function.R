# Purpose: Estimate the mean cumulative function. 

#' Tabulate Events
#' 
#' Tabulate the number at risk and the number of events at each unique
#' observed event time. 
#' 
#' @param time Event time.
#' @param status Status, coded as 1 for the recurrent event, 0 otherwise.
#' @param idx Subject index. 
#' @return Data.frame with these columns:
#' \itemize{
#'   \item `event_times` distinct recurrent event times.
#'   \item `nar` number at risk at the event time. 
#'   \item `events` number of events that occurred. 
#' }

TabulateEvents <- function(time, status, idx){
  
  # Event times.
  event_times <- sort(time[status == 1])
  event_times_tab <- table(event_times)
  unique_event_times <- unique(event_times)
  n_event_times <- length(unique_event_times)
  
  # All censoring times.
  cens_times <- sort(time[status == 0])
  
  # Number initially at risk. 
  init_nar <- length(unique(idx))
  
  # Counts the number of distinct subjects remaining in the data set
  # as of each event time. 
  nar <- sapply(unique_event_times, FUN = function(time){
    return(init_nar - sum(cens_times < time))
  })
  
  # output.
  out <- data.frame(
    'event_times' = unique_event_times, 
    'nar' = nar, 
    'events' = as.numeric(event_times_tab)
  )
  return(out)
}


#' Calculate Mean Cumulative Function
#' 
#' Tabulates the mean cumulative function. See equation 2.1 of 
#' <doi:10.1111/j.0006-341X.2000.00554.x>.
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#' @param idx Unique subject index. 
#' @import survival
#' @export 
#' @return Data.frame with these columns:
#' \itemize{
#'   \item `event_times` distinct recurrent event times.
#'   \item `nar` number at risk at the event time. 
#'   \item `events` number of events that occurred. 
#'   \item `surv`, survival probability.
#'   \item `event_rate`, instantaneous event rate.
#'   \item `mcf`, mean cumulative function. 
#' }

CalculateMCF <- function(
  time,
  status,
  idx
) {
  
  # Form data.
  data <- data.frame(time, status, idx)
  
  # Kaplan-Meier estimator of overall survival. 
  
  # Recode death as 1 and all other events as zero. Keep only the last
  # event for each subject.
  
  # If no deaths are present, the the KM estimator is 1 throughout the 
  # follow-up period. 
  n_death <- sum(data$status == 2)
  if(n_death == 0) {
    shat <- function(x) {1}
  } else {
    
    # Otherwise, estimate the survival function (using deaths only). 
    km_data <- data
    km_data$status <- 1 * (data$status == 2)
    km_data <- km_data[order(km_data$idx, km_data$time, decreasing = c(FALSE, TRUE)), ]
    km_data <- km_data[!duplicated(km_data$idx), ]
    km <- survfit(Surv(time, status) ~ 1, data = km_data)
    shat <- stepfun(x = km$time, y = c(1, km$surv), right = TRUE)
  }

  
  # Nelson-Aalen estimator of event rate.
  rate_data <- data[data$status != 2, ]
  event_tab <- TabulateEvents(rate_data$time, rate_data$status, rate_data$idx)
  event_tab$event_rate <- event_tab$events / event_tab$nar
  
  # Add survival probability.
  event_tab$surv <- shat(event_tab$event_times)
  
  # Estimate marginal cumulative function.
  event_tab$mcf <- cumsum(event_tab$surv * event_tab$event_rate)
  return(event_tab)
}

