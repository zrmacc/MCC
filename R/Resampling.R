#' Grouped Bootstrap.
#'
#' Bootstrap observations grouped by an index. The index `idx` is relabeled
#' in the returned data.frame so that each resampled group has a distinct value.  
#'
#' @param data Data.frame containing `idx`.
#' @return Bootstrapped data.frame.

GroupBoot <- function (data) {
  ids <- sort(unique(data$idx))
  n <- length(ids)

  # Split data by ID.
  split_data <- split(data, data$idx)

  # Sample data.
  key <- sample(x = n, size = n, replace = TRUE)
  out <- split_data[key]

  # Relabel ID.
  Relab <- function(i) {
    sub <- out[[i]]
    sub$idx <- ids[i]
    return(sub)
  }
  out <- lapply(seq_len(n), Relab)

  # Output data.frame.
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  return(out)
}


# -----------------------------------------------------------------------------

#' Stratified, Grouped Bootstrap.
#' 
#' Perform grouped bootstrap within levels of a stratification factor.
#' 
#' @param data Data.frame containing `idx` and `strata`.
#' @return Bootstrapped data.frame.

StratGroupBoot <- function(data) {
  out <- lapply(split(x = data, f = data$strata), GroupBoot)
  out <- do.call(rbind, out)
  return(out)
}


# -----------------------------------------------------------------------------

#' Permute Treatment Assignments.
#'
#' @param data Data.frame containing 'arm' and 'idx'.
#' @importFrom stats rbinom
#' @return Bootstrapped data.frame.

PermData <- function(data) {
  
  # Observed assignments.
  obs_assignments <- unique(data[, c('idx', 'arm')])
  n <- nrow(obs_assignments)
  
  # Permute assignments.
  perm_assignments <- obs_assignments
  perm_assignments$arm <- obs_assignments$arm[sample(n, n, FALSE)]
  perm_trt_arm <- perm_assignments$idx[perm_assignments$arm == 1]
  
  # Format output.
  data$arm <- 0
  data$arm[data$idx %in% perm_trt_arm] <- 1
  return(data)
}
