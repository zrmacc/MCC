#' Generate Bootstrap Data.frame.
#'
#' Bootstrap observations grouped by an index. The index `idx` is relabeled
#' in the returned data.frame so that each resampled group has a distinct value.  
#'
#' @param data Data.frame.
#' @param idx_offset Index offset.
#' @return Bootstrapped data.frame.

BootData <- function(data, idx_offset = 0) {
  ids <- sort(unique(data$idx))
  n <- length(ids)

  # Split data by ID.
  split_data <- split(data, data$idx)

  # Sample data.
  key <- sample(x = n, size = n, replace = TRUE)
  out <- split_data[key]

  # Relabel ID.
  out <- lapply(1:n, function(i) {
    sub <- out[[i]]
    sub$idx <- i + idx_offset 
    return(sub)
  })

  # Output data.frame.
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  return(out)
}

#' Generate Permuted Data.frame.
#'
#' Permute observations grouped by an index by flipping the treatment arm
#' `arm` of an observations with probability 0.5. Arm is assumed to be labeled
#' 1 for treatment, 0 for reference. 
#'
#' @param data Data.frame.
#' @importFrom stats rbinom
#' @return Bootstrapped data.frame.

PermData <- function(data) {
  obs_per_subj <- table(data$idx)
  subj <- length(obs_per_subj)
  
  # Randomize treatment assignment.
  flip <- rbinom(n = subj, size = 1, prob = 0.5)
  flip <- rep(x = flip, times = obs_per_subj)
  
  data$arm <- data$arm * (1 - flip) + (1 - data$arm) * flip
  return(data)
}
