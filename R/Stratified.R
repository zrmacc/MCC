# -----------------------------------------------------------------------------
# AUCs.
# -----------------------------------------------------------------------------

#' Calculate Marginal Area
#' 
#' @param areas Estimated statistics.
#' @param ses Standard errors.
#' @param weights Weights.
#' @param n Sample sizes.
#' @return Data.frame containing:
#' \itemize{
#'   \item The marginal area.
#'   \item The standard error of the area.
#' }

MargAUC <- function(
  areas,
  ses,
  weights,
  n
) {
  out <- data.frame(
    "N" = sum(n),
    "Area" = sum(weights * areas),
    "SE" = sqrt(sum(weights^2 * ses^2))
  )
  return(out)
}


#' Contrast Summary Statistics
#' 
#' Finds the difference and ratio of summary statistics, comparing two arms.
#' 
#' @param area1 Statistic for arm 1.
#' @param area0 Statistic for arm 0.
#' @param se1 Standard error for arm 1.
#' @param se0 Standard error for arm 0.
#' @param alpha Type I error.
#' @importFrom stats pnorm qnorm
#' @return Data.frame containing:
#' \itemize{
#'   \item 'Contrast' and estimate 'Est'.
#'   \item Lower 'L' and upper 'U' confidence bounds.
#'   \item 'P' value.
#' }

ContrastAreas <- function(
  area1,
  area0,
  se1,
  se0,
  alpha 
) {
  crit <- qnorm(p = 1 - alpha / 2)
  
  # Difference.
  delta <- area1 - area0
  se_diff <- sqrt(se1^2 + se0^2)
  delta_lower <- delta - crit * se_diff
  delta_upper <- delta + crit * se_diff
  delta_p <- 2 * pnorm(q = abs(delta) / se_diff, lower.tail = FALSE)
  
  # Ratio.
  rho <- area1 / area0 
  se_rho_log <- sqrt(se1^2 / area1^2 + se0^2 / area0^2)
  rho_lower <- rho * exp(- crit * se_rho_log)
  rho_upper <- rho * exp(+ crit * se_rho_log)
  rho_p <- 2 * pnorm(q = abs(log(rho)) / se_rho_log, lower.tail = FALSE)
  
  # Output.
  out <- data.frame(
    "Contrast" = c("A1-A0", "A1/A0"),
    "Observed" = c(delta, rho),
    "SE" = c(se_diff, rho * se_rho_log),
    "Lower" = c(delta_lower, rho_lower),
    "Upper" = c(delta_upper, rho_upper),
    "P" = c(delta_p, rho_p)
  )
  return(out)
}


#' Calculate Test Statistics
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param strata Optional stratification factor.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param return_areas Return the AUCs?
#' @import dplyr
#' @export 
#' @return If `return_areas`, list containing:
#' \itemize{
#'   \item 'avg_mcf', average MCF across strata.
#'   \item 'contrasts', including the difference and ratio of areas.
#'   \item 'marg_areas', the marginal areas for each arm.
#'   \item 'stratum_areas', the per-stratum and per-arm areas.
#'   \item 'weights', the weights of each stratum.
#' }
#'  Else, only 'contrasts' is returned. 

AUC.Stats.Strat <- function(
  time,
  status,
  arm,
  idx,
  strata,
  tau,
  alpha,
  return_areas = FALSE
) {
  
  # Form data.frame.
  data <- data.frame(
    "time" = time,
    "status" = status,
    "arm" = arm,
    "idx" = idx,
    "strata" = strata
  )
  
  # Stratum sizes.
  stratum_sizes <- data %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise("n" = length(unique(idx)), .groups = "drop") %>%
    as.data.frame
  stratum_sizes$stratum_weight <- stratum_sizes$n / sum(stratum_sizes$n)
  
  # Stratum areas.
  stratum_areas <- data %>%
    dplyr::group_by(arm, strata) %>%
    dplyr::summarise(
      AUC.Area(time, status, idx, tau = tau),
      .groups = "drop"
    ) %>% 
    dplyr::inner_join(
      stratum_sizes[, c("strata", "stratum_weight")], 
      by = "strata"
    ) %>%
    as.data.frame
  
  # Marginal areas.
  area <- se_area <- stratum_weight <- NULL
  marg_areas <- stratum_areas %>%
    dplyr::group_by(arm) %>%
    dplyr::summarise(
      MargAUC(areas = area, ses = se_area, weights = stratum_weight, n = n),
      .groups = "drop"
    ) %>%
    as.data.frame
  marg_areas$Tau <- tau
  colnames(marg_areas)[1] <- "Arm"
  col_order <- c(
    "Arm",
    "N",
    "Tau",
    "Area",
    "SE"
  )
  marg_areas <- marg_areas[, col_order]
  
  # Difference and ratio.
  contrasts <- ContrastAreas(
    area1 = marg_areas$Area[marg_areas$Arm == 1],
    se1 = marg_areas$SE[marg_areas$Arm == 1],
    area0 = marg_areas$Area[marg_areas$Arm == 0],
    se0 = marg_areas$SE[marg_areas$Arm == 0],
    alpha = alpha
  )
  
  # Output
  if(return_areas){
    
    # Average curves.
    mcf1 <- data %>%
      dplyr::filter(arm == 1) %>%
      dplyr::group_by(strata) %>%
      dplyr::summarise(
        CalcMCF(time, status, idx),
        .groups = "keep"
      ) %>%
      dplyr::group_split()
    avg_mcf1 <- AvgMCF(mcf1, weights = stratum_sizes$stratum_weight)
    avg_mcf1$arm <- 1
    
    mcf0 <- data %>%
      dplyr::filter(arm == 0) %>%
      dplyr::group_by(strata) %>%
      dplyr::summarise(
        CalcMCF(time, status, idx),
        .groups = "keep"
      ) %>%
      dplyr::group_split()
    avg_mcf0 <- AvgMCF(mcf1, weights = stratum_sizes$stratum_weight)
    avg_mcf0$arm <- 0
    
    avg_mcf <- rbind(avg_mcf1, avg_mcf0)
    
    # Outputs.
    out <- list(
      'avg_mcf' = avg_mcf,
      'contrasts' = contrasts,
      'marg_areas' = marg_areas,
      'stratum_areas' = stratum_areas
    )
  } else {
    out <- contrasts
  }
  return(out)
}


# -----------------------------------------------------------------------------
# Bootstrap/permutation
# -----------------------------------------------------------------------------

#' Bootstrap Inference
#'
#' Constructs bootstrap confidence intervals.
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param strata Optional stratification factor.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @param seed Simulation seed.
#' @return Data.frame containing:
#' \itemize{
#'   \item Bootstrap difference 'boot_diff' and ratio 'boot_ratio' of areas.
#'   \item An indicator that the bootstrap difference was of the opposite
#'     sign, 'is_diff_sign'.
#' }

Boot.Sim.Strat <- function(
  time, 
  status,
  arm,
  idx,
  strata,
  obs_stats,
  tau,
  alpha,
  reps,
  seed
) {
  set.seed(seed)
  data <- data.frame(
    time = time,
    status = status,
    arm = arm,
    idx = idx,
    strata = strata
  )
  data1 <- subset(x = data, arm == 1)
  data0 <- subset(x = data, arm == 0)
  n0 <- nrow(data0)

  # Bootstrap function.
  loop <- function(b) {

    # Bootstrap data sets.
    boot0 <- StratGroupBoot(data0, idx_offset = 0)
    boot1 <- StratGroupBoot(data1, idx_offset = n0)
    boot <- rbind(boot0, boot1)

    # Bootstrap statistics.
    boot_stats <- AUC.Stats.Strat(
      time = boot$time,
      status = boot$status,
      arm = boot$arm,
      idx = boot$idx,
      strata = boot$strata,
      tau = tau,
      alpha = alpha
    )
    names(boot_stats) <- paste0("Boot_", names(boot_stats))

    # Bootstrap p-value indicators.
    # Indicator is 1 if the sign of the difference in areas is opposite that observed.
    is_diff_sign <- sign(boot_stats$Boot_Observed[1]) != sign(obs_stats$Observed[1])

    # Results
    out <- c(
      boot_stats$Boot_Observed,
      is_diff_sign
    )
    return(out)
  }

  sim <- lapply(seq(1:reps), loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c("boot_diff", "boot_ratio", "is_diff_sign")
  return(sim)
}


#' Bootstrap Confidence Intervals
#'
#' @param sim Bootstrap simulation results, as generated by \code{\link{Boot.Sim.Strat}}.
#' @param obs_stats Observed contrasts.
#' @param alpha Type I error.
#' @importFrom stats sd
#' @return Data.frame containing the equi-tailed and highest-density bootstrap
#'   confidence intervals.

Boot.CIs.Strat <- function(
  sim,
  obs_stats,
  alpha
) {

  # Equi-tailed CI for difference.
  eti_diff <- HighDensCI(
    x = sim$boot_diff,
    min_tail_prob = alpha / 2,
    intervals = 0
  )

  # Equi-tailed CI for ratio.
  eti_ratio <- HighDensCI(
    x = log(sim$boot_ratio),
    min_tail_prob = alpha / 2,
    intervals = 0
  )
  eti_ratio[1:2] <- exp(eti_ratio[1:2])

  # HDI for difference.
  reps <- nrow(sim)
  hdi_diff <- HighDensCI(
    x = sim$boot_diff,
    alpha = alpha,
    min_tail_prob = 1 / reps,
    intervals = 1e3
  )

  # HDI for ratio.
  hdi_ratio <- HighDensCI(
    x = log(sim$boot_ratio),
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
  se_diff <- sd(sim$boot_diff)
  sim$log_boot_ratio <- log(sim$boot_ratio)
  se_ratio <- exp(mean(sim$log_boot_ratio)) * sd(sim$log_boot_ratio)
  rownames(cis) <- NULL
  cis$Method <- "Bootstrap"
  cis$Type <- rep(c("Equitailed", "Highest-density"), times = 2)
  cis$Contrast <- rep(c("A1-A0", "A1/A0"), each = 2)
  cis$Observed <- rep(obs_stats$Observed, each = 2)
  cis$SE <- rep(c(se_diff, se_ratio), each = 2)
  col_order <- c(
    "Method", 
    "Type", 
    "Contrast", 
    "Observed",
    "SE",
    "Lower", 
    "Upper", 
    "Alpha_Lower", 
    "Alpha_Upper"
  )
  cis <- cis[, col_order]
  return(cis)
}


# -----------------------------------------------------------------------------

#' Permutation Inference
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param strata Optional stratification factor.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @param seed Simulation seed.
#' @return Data.frame containing:
#' \itemize{
#'   \item Permutation difference 'perm_diff' and ratio 'perm_ratio' of areas.
#'   \item Indicators that the permutation difference or ratio was as or more extreme
#'     than the observed difference or ratio.
#' }

CompAUCs.Perm.Strat <- function(
  time,
  status,
  arm,
  idx,
  strata,
  obs_stats,
  tau,
  alpha,
  reps,
  seed
) {
  set.seed(seed)
  data <- data.frame(
    time = time,
    status = status,
    arm = arm,
    idx = idx,
    strata = strata
  )

  # Permutation function.
  loop <- function(b) {

    # Permute data.
    perm <- PermData(data)

    # Permutation statistics
    perm_stats <- AUC.Stats.Strat(
      time = perm$time,
      status = perm$status,
      arm = perm$arm,
      idx = perm$idx,
      strata = perm$strata,
      tau = tau,
      alpha = alpha
    )
    names(perm_stats) <- paste0("Perm_", names(perm_stats))

    # Permutation indicators.
    perm_diff <- perm_stats$Perm_Observed[1]
    obs_diff <- obs_stats$Observed[1]
    perm_ratio <- perm_stats$Perm_Observed[2]
    obs_ratio <- obs_stats$Observed[2]

    is_diff_sign <- (sign(perm_diff) != sign(obs_diff))
    is_diff_more_extreme <- abs(perm_diff) >= abs(obs_diff)
    is_ratio_more_extreme <- abs(log(perm_ratio)) >= abs(log(obs_ratio))
    perm_diff_1sided <- is_diff_sign * is_diff_more_extreme
    perm_ratio_1sided <- is_diff_sign * is_ratio_more_extreme

    # Results
    out <- c(
      perm_stats$Perm_Observed,
      "perm_diff_1sided" = perm_diff_1sided,
      "perm_ratio_1sided" = perm_ratio_1sided
    )
    return(out)
  }

  sim <- lapply(seq(1:reps), loop)
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
# Main function
# -----------------------------------------------------------------------------

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
#' @param alpha Alpha level.
#' @param boot Logical, construct bootstrap confidence intervals?
#' @param perm Logical, perform permutation test?
#' @param reps Replicates for bootstrap/permutation inference.
#' @param seed Seed for bootstrap/permutation inference.
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
#' # Simulate data set.
#' data <- GenData()
#'
#' aucs <- CompareStratAUCs(
#'   time = data$time,
#'   status = data$status,
#'   arm = data$arm,
#'   idx = data$idx,
#'   strata = data$strata,
#'   tau = 2,
#'   boot = TRUE,
#'   perm = TRUE,
#'   reps = 25,
#'   alpha = 0.05
#' )
#' show(aucs)
#' }

CompareStratAUCs <- function(
  time,
  status,
  arm,
  idx,
  tau,
  strata = NULL,
  alpha = 0.05,
  boot = FALSE,
  perm = FALSE,
  reps = 2000,
  seed = 100
) {

  # Create single stratum if no strata are provided.
  if (is.null(strata)) {
    strata <- rep(1, length(time))
  }

  # Observed test stats.
  obs <- AUC.Stats.Strat(
    time = time,
    status = status,
    arm = arm,
    idx = idx,
    strata = strata,
    tau = tau,
    alpha = alpha,
    return_areas = TRUE
  )
  obs_stats <- obs$contrasts

  # CIs.
  cis <- obs_stats[, c("Contrast", "Observed", "SE", "Lower", "Upper")]
  cis <- data.frame(
    "Method" = "Asymp",
    "Type" = "Equitailed",
    cis,
    "Alpha_Lower" = alpha / 2,
    "Alpha_Upper" = alpha / 2
  )

  # P-values.
  pvals <- obs_stats[, c("Contrast", "Observed", "P")]
  pvals <- data.frame(
    "Method" = "Asymp",
    pvals
  )

  # Simulation replicates.
  sim_reps <- list()

  # -------------------------------------------------------

  # Bootstrap inference.
  if (boot) {

    # Simulate.
    boot_sim <- Boot.Sim.Strat(
      time = time,
      status = status,
      arm = arm,
      idx = idx,
      strata = strata,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps,
      seed = seed
    )
    sim_reps$boot_sim <- boot_sim

    # Confidence intervals.
    boot_cis <- Boot.CIs.Strat(
      sim = boot_sim,
      obs_stats = obs_stats,
      alpha = alpha
    )
    cis <- rbind(
      cis,
      boot_cis
    )
    cis <- cis[order(cis$Contrast), ]

    # P-value.
    boot_p <- min(2 * mean(c(1, boot_sim$is_diff_sign)), 1)
    boot_pvals <- cbind(
      "Method" = "Bootstrap",
      pvals[, c("Contrast", "Observed")],
      "P" = boot_p
    )
    pvals <- rbind(
      pvals,
      boot_pvals
    )
  }

  # -------------------------------------------------------

  # Permutation inference.
  if (perm) {

    # Simulate.
    perm_sim <- CompAUCs.Perm.Strat(
      time = time,
      status = status,
      arm = arm,
      idx = idx,
      strata = strata,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps,
      seed = seed
    )
    sim_reps$perm_sim <- perm_sim

    # Permutation p-values.
    perm_pvals <- apply(perm_sim[, c(3:4)], 2, function(x) {
      return(min(2 * mean(c(1, x)), 1))
    })
    perm_pvals <- data.frame(
      "Method" = "Perm",
      "Contrast" = c("A1-A0", "A1/A0"),
      "Observed" = obs_stats$Observed,
      "P" = perm_pvals
    )
    rownames(perm_pvals) <- NULL
    pvals <- rbind(
      pvals,
      perm_pvals
    )
    pvals <- pvals[order(pvals$Observed), ]
  }

  # -------------------------------------------------------

  # Output
  out <- new(
    Class = "CompStratAUCs",
    StratumAreas = obs$stratum_areas,
    MargAreas = obs$marg_areas,
    CIs = cis,
    MCF = obs$avg_mcf,
    Reps = sim_reps,
    Pvals = pvals
  )

  return(out)
}
