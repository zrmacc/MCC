#' Inference on the Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0).
#' 
#' Two methods of p-value calculation are available. For 'perm', treatment 
#' assignments are permuted on each iteration, and the p-value is the
#' proportion of the *null* statistics that are as or more extreme than
#' the *observed* statistics. For 'boot', the p-value is twice the proportion
#' of bootstrap replicates on which the sign of the difference is areas is 
#' reversed. 
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period. 
#' @param arm Arm, coded as 1 for treatment, 0 for reference. 
#' @param idx Unique subject index. 
#' @param tau Truncation time. 
#' @param strata Optional stratification factor. 
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom stats quantile 
#' @importFrom methods new
#' @export 
#' @return Object of class compAUCs with these slots:
#' \itemize{
#'   \item `@Areas`: The AUC for each arm.
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals. 
#'   \item `@Curves`: Mean cumulative count curve for each arm; averaged across strata 
#'     if present. 
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics. 
#'   \item `@Weights`: Per-stratum weights and AUCs. 
#' }
#' @examples 
#' \donttest{
#' data(mcc_data)
#' set.seed(100)
#' aucs <- CompareAUCs(
#'   time = mcc_data$time,
#'   status = mcc_data$status,
#'   arm = mcc_data$arm,
#'   idx = mcc_data$idx,
#'   tau = 60,
#'   strata = mcc_data$strata,
#'   reps = 100,
#'   alpha = 0.05
#' )
#' show(aucs)
#' }

CompareAUCs <- function(
  time, 
  status, 
  arm, 
  idx, 
  tau, 
  strata = NULL,
  reps = 2000, 
  alpha = 0.05
) { 
  
  # Create single stratum if no strata are provided. 
  if (is.null(strata)){
    strata <- rep(1, length(time))
  }
  
  # Form data.
  data <- data.frame(time, status, arm, idx, strata)
  data$strata <- factor(data$strata)
  
  # Split data. 
  data0 <- data[data$arm == 0, ]
  data1 <- data[data$arm == 1, ]
  n0 <- nrow(data0)
  
  # Observed test stats.
  obs <- AUC.Stats(
    data0 = data0,
    data1 = data1,
    tau = tau, 
    return_areas = TRUE
  ) 
  obs_stats <- obs$stats
  
  # -------------------------------------------------------
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot0 <- StratGroupBoot(data0, idx_offset = 0)
    boot1 <- StratGroupBoot(data1, idx_offset = n0)
    
    # Bootstrap statistics.
    boot_stats <- AUC.Stats(
      data0 = boot0,
      data1 = boot1,
      tau = tau
    )
    names(boot_stats) <- paste0('boot_', names(boot_stats))
    
    # Bootstrap p-value indicators.
    # Indicator is 1 if the sign of the difference in areas is opposite that observed.
    ind_1side_boot <- as.numeric(sign(boot_stats[1]) != sign(obs_stats[1]))
      
    # Permute data.
    perm <- PermData(data)
    perm0 <- perm[perm$arm == 0, ]
    perm1 <- perm[perm$arm == 1, ]
    
    # Permutation statistics
    perm_stats <- AUC.Stats(
      data0 = perm0,
      data1 = perm1,
      tau = tau
    )
    names(perm_stats) <- paste0('perm_', names(perm_stats))
    
    # Permutation indicators.
    sign_diff <- 1 * (sign(perm_stats[1]) != sign(obs_stats[1]))
    ind_2side_perm_diff <- as.numeric(abs(perm_stats[1]) >= abs(obs_stats[1]))
    ind_2side_perm_ratio <- as.numeric(abs(log(perm_stats[2])) >= abs(log(obs_stats[2])))
    ind_1side_perm_diff <- sign_diff * ind_2side_perm_diff
    ind_1side_perm_ratio <- sign_diff * ind_2side_perm_ratio
    
    # Results
    out <- c(
      boot_stats, 
      perm_stats,
      '1side_boot' = ind_1side_boot,
      '1side_boot' = ind_1side_boot,
      '1side_perm_diff' = ind_1side_perm_diff,
      '2side_perm_diff' = ind_2side_perm_diff,
      '1side_perm_ratio' = ind_1side_perm_ratio,
      '2side_perm_ratio' = ind_2side_perm_ratio
    )
    return(out)
  }
  
  sim <- lapply(seq(1:reps), aux)
  sim <- do.call(rbind, sim)
  
  # -------------------------------------------------------
  
  # Equi-tailed CI for difference.
  eti_diff <- HighDensCI(
    x = sim[, 1], 
    min_tail_prob = alpha / 2,
    intervals = 0
  )
  
  # Equi-tailed CI for ratio.
  eti_ratio <- HighDensCI(
    x = log(sim[, 2]),
    min_tail_prob = alpha / 2,
    intervals = 0
  )
  eti_ratio[1:2] <- exp(eti_ratio[1:2])
  
  # HDI for difference.
  hdi_diff <- HighDensCI(
    x = sim[, 1],
    alpha = alpha,
    min_tail_prob = 1 / reps,
    intervals = 1e3
  )
  
  # HDI for ratio.
  hdi_ratio <- HighDensCI(
    x = log(sim[, 2]),
    alpha = alpha,
    min_tail_prob = 1 / reps,
    intervals = 1e3
  )
  hdi_ratio[1:2] <- exp(hdi_ratio[1:2])
  
  # Format confidence intervals. 
  cis <- data.frame(
    rbind(
      eti_diff,
      hdi_diff,
      eti_ratio,
      hdi_ratio
    )
  )
  rownames(cis) <- NULL
  colnames(cis) <- c('Lower', 'Upper', 'Lower_alpha', 'Upper_alpha')
  cis$Method <- rep(c('Equi-tail', 'Highest-density'), times = 2)
  cis$Contrast <- rep(c('A1-A0', 'A1/A0'), each = 2)
  cis$Observed <- rep(obs_stats, each = 2)
  cis <- cis[, c(5:7, 1:4)]
  
  # -------------------------------------------------------
  
  # P-values.
  pvals <- data.frame(
    'Method' = c(rep('Boot', times = 2), rep('Perm', times = 4)),
    'Sides' = c(rep(1, times = 2), rep(c(1, 2), times = 2)),
    'Contrast' = c(c('A1-A0', 'A1/A0'), rep(c('A1-A0', 'A1/A0'), each = 2)),
    'P' = apply(sim[, 5:10], 2, mean)
  )
  
  # -------------------------------------------------------
  
  # Output
  out <- new(
    Class = 'compAUCs',
    Areas = obs$areas,
    CIs = cis,
    Curves = obs$mcf,
    Reps = sim[, 1:4],
    Pvals = pvals,
    Weights = obs$weights
  )
  
  return(out)
}

