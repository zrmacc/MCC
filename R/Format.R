# Purpose: Format data.

#' Format Data for a Single Subject
#' 
#' @param df Data.frame for a single subject.
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @return Data.frame with an added censoring event.
#' @noRd

FormatSubj <- function(df, cens_after_last) {
  obs_end <- sum(df$obs_end)
  out <- df
  
  # Check for multiple observation terminating events.
  if (obs_end > 1) {
    stop(paste0("Subject ", unique(df$idx), " has multiple observation terminating events."))
  }
  
  # Add censoring if no observation terminating event is present.
  if (obs_end == 0 & cens_after_last) {
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
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @return Formatted data.frame.
#' @importFrom dplyr "%>%"
#' @export 

FormatData <- function(
  idx,
  time,
  status,
  arm, 
  covars,
  strata,
  cens_after_last
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
  data <- data %>%
    dplyr::mutate(obs_end = 1 * (status != 1)) %>%
    dplyr::arrange(idx, time, obs_end)
  
  # Format data.
  split_data <- split(x = data, f = data$idx)
  format_data <- lapply(split_data, function(x) {
    FormatSubj(x, cens_after_last = cens_after_last)
  })
  final_data <- do.call(rbind, format_data)
  rownames(final_data) <- NULL
  
  # Ensure each subject has exactly 1 terminal event.
  check <- final_data %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      n_obs_end = sum(obs_end)
    )
  
  if (any(check$n_obs_end != 1)) {
    warning("Patience without censoring times were found.")
  }
  
  final_data$obs_end <- NULL
  return(final_data)
}