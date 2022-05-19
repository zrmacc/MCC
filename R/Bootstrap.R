# Purpose: Functions for bootstrap simulation.

# -----------------------------------------------------------------------------
# Bootstrap/permutation
# -----------------------------------------------------------------------------

#' Bootstrap Inference for Stratified Estimator
#'
#' Constructs bootstrap confidence intervals.
#'
#' @param data Data.frame containing: idx, time, status, arm, strata.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return Data.frame containing:
#' \itemize{
#'   \item Bootstrap difference 'boot_diff' and ratio 'boot_ratio' of areas.
#'   \item An indicator that the bootstrap difference was of the opposite
#'     sign, 'is_diff_sign'.
#' }

BootSimStrat <- function(
  data,
  obs_stats,
  tau,
  alpha = 0.05,
  reps = 2000
) {
  
  # Partition data.
  arm <- NULL
  data1 <- subset(x = data, arm == 1)
  data0 <- subset(x = data, arm == 0)
  
  # Bootstrap function.
  Loop <- function(b) {
    
    # Bootstrap data sets.
    boot0 <- StratGroupBoot(data0)
    boot1 <- StratGroupBoot(data1)
    boot <- rbind(boot0, boot1)
    
    # Bootstrap statistics.
    boot_stats <- CalcStratAUC(
      data = boot,
      tau = tau,
      alpha = alpha,
      return_areas = FALSE
    )
    names(boot_stats) <- paste0("boot_", names(boot_stats))
    
    # Bootstrap p-value indicators.
    # Indicator is 1 if the sign of the difference in areas is opposite that observed.
    is_diff_sign <- sign(boot_stats$boot_observed[1]) != sign(obs_stats$observed[1])
    
    # Results
    out <- c(
      boot_stats$boot_observed,
      is_diff_sign
    )
    return(out)
  }
  
  sim <- lapply(seq_len(reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c("boot_diff", "boot_ratio", "is_diff_sign")
  return(sim)
}


# -----------------------------------------------------------------------------


#' Bootstrap Inference for Augmentation Estimator
#'
#' Constructs bootstrap confidence intervals.
#'
#' @param data Data.frame containing: idx, time, status, arm, strata.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return Data.frame containing:
#' \itemize{
#'   \item Bootstrap difference 'boot_diff' and ratio 'boot_ratio' of areas.
#'   \item An indicator that the bootstrap difference was of the opposite
#'     sign, '1side_boot_diff'.
#' }

BootSimAug <- function(
  data,
  obs_stats,
  tau,
  alpha = 0.05,
  reps = 2000
) {
  
  # Partition data.
  arm <- NULL
  data1 <- subset(x = data, arm == 1)
  data0 <- subset(x = data, arm == 0)
  
  # Bootstrap function.
  Loop <- function(b) {
    
    # Bootstrap data sets.
    boot0 <- GroupBoot(data0)
    boot1 <- GroupBoot(data1)
    boot <- rbind(boot1, boot0)
    
    # Bootstrap statistics.
    boot_stats <- CalcAugAUC(
      data = boot,
      tau = tau,
      alpha = alpha
    )
    names(boot_stats) <- paste0("boot_", names(boot_stats))
    
    # Bootstrap p-value indicators.
    # Indicator is 1 if the sign of the difference in areas is opposite that observed.
    is_diff_sign <- sign(boot_stats$boot_observed[1]) != sign(obs_stats$observed[1])
    
    # Results
    out <- c(
      boot_stats$boot_observed,
      is_diff_sign
    )
    return(out)
  }
  
  sim <- lapply(seq_len(reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c("boot_diff", "is_diff_sign")
  return(sim)
}