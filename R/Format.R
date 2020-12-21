# Purpose: Format data.

#' Format Data for a Single Subject
#' 
#' @param df Data.frame for a single subject.
#' @return Data.frame with an added censoring event.

FormatSubj <- function(df) {
  obs_end <- sum(df$obs_end)
  out <- df
  
  # Check for multiple observation terminating events.
  if (obs_end > 1) {
    stop(paste0("Subject ", unique(df$idx), " has multiple observation terminating events."))
  }
  
  # Add censoring if no observation terminating event is present.
  if (obs_end == 0) {
    n_row <- nrow(df)
    last_row <- df[n_row, ]
    last_row$status <- 0
    last_row$obs_end <- 1
    out <- rbind(out, last_row)
  }
  
  return(out)
}

#' Format Data
#' 
#' @param idx Unique subject index.
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
#' @param strata Optional stratification factor.
#' @return Formatted data.frame.

FormatData <- function(
  idx,
  time,
  status,
  arm, 
  covars,
  strata
) {
  
  data <- data.frame(idx, time, status, arm)
  
  if (is.null(covars) & is.null(strata)) {
    data$strata <- 1
  } else if (!is.null(covars) & is.null(strata)) {
    data <- cbind(data, covars)
  } else if(is.null(covars) & !is.null(strata)) {
    data$strata <- strata
  }
  
  # Indicator of an observation-terminating event.
  data$obs_end <- 1 * (data$status != 1)
  
  # Sort data.
  data <- data[order(data$idx, data$time, data$obs_end), ]
  
  # Format data.
  split_data <- split(x = data, f = data$idx)
  format_data <- lapply(split_data, FormatSubj)
  final_data <- do.call(rbind, format_data)
  rownames(final_data) <- NULL
  
  # Ensure each subject has exactly 1 terminal event.
  check <- tapply(final_data$obs_end, final_data$idx, sum)
  if (unique(check) != 1){
    stop("Error formatting the data. Please ensure each subject has a single observation terminating event.")
  }
  
  final_data$obs_end <- NULL
  return(final_data)
}