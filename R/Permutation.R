# Purpose: Functions for permutation inference.
# Updated: 2022-05-19

# -----------------------------------------------------------------------------
# P-value calculation.
# -----------------------------------------------------------------------------

#' Calculate P
#' 
#' Calculates p-value from 1-sided rejection indicator.
#' @param p Vector of 1/0 rejection indicators.
#' @return Numeric p-value.
#' @noRd

CalcP <- function(p) {
  p1 <- c(1, p)
  mu <- mean(p1)
  if (mu <= 0.5) {
    out <- 2 * mu
  } else {
    out <- 2 * mean(c(1, 1 - p))
  }
  out <- min(out, 1)
  return(out)
}

# -----------------------------------------------------------------------------
# Stratified Estimator.
# -----------------------------------------------------------------------------

#' Permutation Inference for Stratified Estimator
#'
#' @param data Data.frame containing: {arm, idx, status, strata time}.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return Data.frame containing:
#' \itemize{
#'   \item Permutation difference 'perm_diff' and ratio 'perm_ratio' of areas.
#'   \item Indicators that the permutation difference or ratio was as or more extreme
#'     than the observed difference or ratio.
#' }

PermSimStrat <- function(
  data,
  obs_stats,
  tau,
  alpha = 0.05,
  reps = 2000
) {
  
  # Permutation function.
  Loop <- function(b) {
    
    # Permute data.
    perm <- PermData(data)
    
    # Permutation statistics
    perm_stats <- CalcStratAUC(
      data = perm,
      tau = tau,
      alpha = alpha
    )
    names(perm_stats) <- paste0("perm_", names(perm_stats))
    
    # Permutation indicators.
    perm_diff <- perm_stats$perm_observed[1]
    obs_diff <- obs_stats$observed[1]
    perm_ratio <- perm_stats$perm_observed[2]
    obs_ratio <- obs_stats$observed[2]
    
    is_diff_sign <- (sign(perm_diff) != sign(obs_diff))
    is_diff_more_extreme <- abs(perm_diff) >= abs(obs_diff)
    is_ratio_more_extreme <- abs(log(perm_ratio)) >= abs(log(obs_ratio))
    perm_diff_1sided <- is_diff_sign * is_diff_more_extreme
    perm_ratio_1sided <- is_diff_sign * is_ratio_more_extreme
    
    # Results
    out <- c(
      perm_stats$perm_observed,
      perm_diff_1sided = perm_diff_1sided,
      perm_ratio_1sided = perm_ratio_1sided
    )
    return(out)
  }
  
  sim <- lapply(seq(1:reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c(
    "perm_diff",
    "perm_ratio",
    "perm_diff_1sided",
    "perm_ratio_1sided"
  )
  return(sim)
}


# -----------------------------------------------------------------------------
# Augmentation Estimator.
# -----------------------------------------------------------------------------

#' Permutation Inference for Augmentation Estimator
#'
#' @param data Data.frame containing: idx, time, status, arm, covars.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return Data.frame containing:
#' \itemize{
#'   \item Permutation difference 'perm_diff' and ratio 'perm_ratio' of areas.
#'   \item Indicators that the permutation difference or ratio was as or more extreme
#'     than the observed difference or ratio.
#' }

PermSimAug <- function(
  data,
  obs_stats,
  tau,
  alpha = 0.05,
  reps = 2000
) {
  
  # Permutation function.
  Loop <- function(b) {
    
    # Permute data.
    perm <- PermData(data)
    
    # Permutation statistics
    perm_stats <- CalcAugAUC(
      data = perm,
      tau = tau,
      alpha = alpha
    )
    names(perm_stats) <- paste0("perm_", names(perm_stats))
    
    # Permutation indicators.
    perm_diff <- perm_stats$perm_observed[1]
    obs_diff <- obs_stats$observed[1]
    
    is_diff_sign <- sign(perm_diff) != sign(obs_diff)
    is_more_extreme <- abs(perm_diff) >= abs(obs_diff)
    
    # Results
    out <- c(
      perm_stats$perm_observed,
      "perm_1sided" = is_diff_sign * is_more_extreme
    )
    return(out)
  }
  
  sim <- lapply(seq_len(reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c("perm_diff", "perm_1sided")
  return(sim)
}
