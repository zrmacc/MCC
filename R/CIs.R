# -----------------------------------------------------------------------------

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
    "L" = as.numeric(lower[key]),
    "U" = as.numeric(upper[key]),
    "alpha_L" = lower_probs[key],
    "alpha_U" = 1 - upper_probs[key]
  )
  return(out)
}