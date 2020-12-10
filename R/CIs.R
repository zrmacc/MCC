#' Highest Density Confidence Interval
#' 
#' @param x Bootstrap replicates.
#' @param alpha Alpha level.
#' @param min_tail_prob Minimum alpha per tail. 
#' @param intervals Intervals to consider. 

HighDensCI <- function(
  x,
  alpha = 0.05,
  min_tail_prob = 0.005,
  intervals = 1e2
) {
  
  # Evaluation grid.
  lower_probs <- seq(
    from = min_tail_prob,
    to = alpha - min_tail_prob,
    length.out = intervals + 1
  )
  upper_probs <- rev(1 - lower_probs)
  
  # Bounds.
  lower <- quantile(x = x, probs = lower_probs, na.rm = TRUE)
  upper <- quantile(x = x, probs = upper_probs, na.rm = TRUE)
  
  # Minimum length.
  len <- abs(upper - lower)
  
  # Final interval.
  key <- which.min(len)
  
  # Output.
  out <- c(
    "Lower" = as.numeric(lower[key]),
    "Upper" = as.numeric(upper[key]),
    "Alpha_Lower" = lower_probs[key],
    "Alpha_Upper" = 1 - upper_probs[key]
  )
  return(out)
}