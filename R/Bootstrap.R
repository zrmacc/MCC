#' Generate Bootstrap Data.frame.
#'
#' Bootstrap observations grouped by an index. The index `idx` is relabeled
#' in the returned data.frame so that resampled group has a distinct value.  
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
  for (i in 1:n) {
    out[[i]]$idx <- i + idx_offset
  }

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

PermuteData <- function(data) {
  ids <- sort(unique(data$idx))
  n <- length(ids)

  # Split data by ID.
  split_data <- split(data, data$idx)

  # Relabel ID.
  for (i in 1:n) {
    flip <- rbinom(n = 1, size = 1, prob = 0.5)
    if (flip == 1) {
      split_data[[i]]$arm <- 1 - split_data[[i]]$arm
    }
  }

  # Output data.frame.
  out <- do.call(rbind, split_data)
  rownames(out) <- NULL
  return(out)
}
