#' Grouped Bootstrap.
#'
#' Bootstrap observations grouped by an index. The index `idx` is relabeled
#' in the returned data.frame so that each resampled group has a distinct value.  
#'
#' @param data Data.frame containing `idx`.
#' @param idx_offset Index offset.
#' @return Bootstrapped data.frame.

GroupBoot <- function (data, idx_offset = 0) {
  ids <- sort(unique(data$idx))
  n <- length(ids)

  # Split data by ID.
  split_data <- split(data, data$idx)

  # Sample data.
  key <- sample(x = n, size = n, replace = TRUE)
  out <- split_data[key]

  # Relabel ID.
  out <- lapply(1:n, function (i) {
    sub <- out[[i]]
    sub$idx <- i + idx_offset 
    return(sub)
  })

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
#' @param idx_offset Index offset.
#' @return Bootstrapped data.frame.

StratGroupBoot <- function(data, idx_offset = 0) {
  
  # Partition by strata.
  data_strata <- split(data, data$strata, drop = TRUE)
  
  # Stratum sizes.
  sizes <- sapply(
    data_strata,
    function (x) {length(unique(x$idx))}
  )
  offsets <- idx_offset + c(0, sizes[1:(length(sizes) - 1)])
  
  # Bootstrap within strata. 
  out <- lapply(
    seq_along(data_strata),
    function (i) {GroupBoot(data_strata[[i]], offsets[i])}
  )
  
  # Output
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
