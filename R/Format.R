# Purpose: Format input data.
# Updated: 2025-11-08

#' Format Data for a Single Subject
#' 
#' @param df Data.frame for a *single subject*.
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @return Data.frame with an added censoring event.
#' @noRd
FormatSubj <- function(df, cens_after_last = TRUE) {
  obs_end <- sum(df$obs_end)
  out <- df
  
  # Check for multiple observation terminating events.
  if (obs_end > 1) {
    which_idx <- unique(df$orig_idx)
    stop(glue::glue("Subject {which_idx} has multiple observation terminating events."))
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


#' Map Original Index to an Integer
#' 
#' Data is expected to have an integer index. This function creates
#' an internal index in case the original index is provided in another
#' format.
#'
#' @param data Data.Frame.
#' @return List containing the data.frame with idx converted to an 
#'   integer plus the mapping from the original to the new index.
#' @noRd
ConvertIdxToInt <- function(data) {
  idx <- orig_idx <- NULL
  data$orig_idx <- data$idx
  data$idx <- NULL
  
  idx_data <- data %>% 
    dplyr::select(orig_idx) %>% 
    unique()
  
  idx_data$idx <- seq_len(nrow(idx_data))
  out <- merge(x = data, y = idx_data, by = "orig_idx")
  out <- out %>% dplyr::relocate(idx)
  return(out)
}


#' Format Data
#' 
#' Subsets the key columns for analysis, converts the index to an integer,
#' adds placeholders for `strata` and `weights` if not specified, checks
#' to ensure each unique index has 1 and only 1 o
#' 
#' @param data Data.frame.
#' @param arm_name Name of column containing treatment arm. Must be coded as 1
#'   for treatment, 0 for reference.
#' @param cens_after_last Should subjects who lack an explicit censoring time
#'   be censored after their last observed event? 
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
#' @param idx_name Name of column containing a unique subject index.
#' @param keep_orig_idx Keep the original index? 
#' @param status_name Name of column containing the status. Must be coded as 0 for
#'   censoring, 1 for event, 2 for death. Each subject should have an 
#'   observation-terminating event, either censoring or death.
#' @param strata Optional stratification factor. Should not be provided if a
#'   covariate matrix is provided.
#' @param time_name Name of column containing the observation time.
#' @param weights Optional column of weights.
#' @return Formatted data.frame.
#' @export 
FormatData <- function(
  data,
  arm_name = "arm",
  cens_after_last = TRUE,
  covars = NULL,
  idx_name = "idx",
  keep_orig_idx = FALSE,
  status_name = "status",
  strata = NULL,
  time_name = "time",
  weights = NULL
) {
  
  # Rename columns as necessary.
  if (!is.null(arm_name)) {
    # Two-sample case.
    arm <- idx <- status <- time <- NULL
    key_cols <- c(arm_name, idx_name, status_name, time_name)
    data <- data %>%
      dplyr::select(dplyr::all_of(key_cols)) %>%
      dplyr::rename(
        arm = {{arm_name}},
        idx = {{idx_name}},
        status = {{status_name}},
        time = {{time_name}}
      )
  } else {
    # Single-sample case.
    idx <- status <- time <- NULL
    key_cols <- c(idx_name, status_name, time_name)
    data <- data %>%
      dplyr::select(dplyr::all_of(key_cols)) %>%
      dplyr::rename(
        idx = {{idx_name}},
        status = {{status_name}},
        time = {{time_name}}
      )
  }

  # Ensure index is an integer.
  data <- ConvertIdxToInt(data)
  
  # Add covariates or strata.
  if (is.null(covars) & is.null(strata)) {
    data$strata <- 1
  } else if (!is.null(covars) & is.null(strata)) {
    data <- cbind(data, covars)
  } else if(is.null(covars) & !is.null(strata)) {
    data$strata <- strata
  }
  
  # Add jump weights.
  if (is.null(weights)) {
    data$weights <- 1
  } else {
    data$weights <- weights
  }
  
  # Crate an indicator of an observation-terminating event.
  # Sort by index > time > terminal indicator.
  data <- data %>%
    dplyr::mutate(obs_end = 1 * (status != 1)) %>%
    dplyr::arrange(idx, time, obs_end)
  
  # Format each subject.
  split_data <- split(x = data, f = data$idx)
  format_data <- lapply(split_data, function(x) {
    FormatSubj(x, cens_after_last = cens_after_last)
  })
  data <- do.call(rbind, format_data)
  rownames(data) <- NULL
  
  # Ensure each subject has exactly 1 terminal event.
  obs_end <- NULL
  check <- data %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      n_obs_end = sum(obs_end)
    )
  if (any(check$n_obs_end != 1)) {
    warning("Patients without observation terminating events were found.")
  }
  data$obs_end <- NULL
  
  # Format.
  if (!keep_orig_idx) {
    data$orig_idx <- NULL
  }
  return(data)
}
