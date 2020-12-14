# Purpose: Augmentation estimator.
# Updated: 2020-12-08

#' AUC Calculation for Augmentation Estimator.
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to remain
#'   at risk throughout the observation period.
#' @param idx Unique subject index.
#' @param subj_covars Per-subject covariate data. Should include `idx` and
#'   covariates only.
#' @param mu Grand mean.
#' @param tau Truncation time. 
#' @param calc_var Calculate analytical variance?
#' @return List containing:
#' \itemize{
#'   \item Tabulated mean cumulative function 'mcf'.
#'   \item Estimated 'area' and its variance 'var_area'.
#'   \item 'sigma', \eqn{\frac{1}{n^2}\sum_{i=1}^{n}(X_{i}-\mu)(X_{i}-\mu)'}.
#'   \item 'gamma', \eqn{\frac{1}{n^2}\sum_{i=1}^{n}(X_{i}-\mu)\xi_{i}}.
#'   \item 'mean_resid', \eqn{\frac{1}{n}\sum_{i=1}^{n}(X_{i}-\mu)}.
#' }

AUC.Area.Aug <- function(
  time, 
  status,
  idx,
  subj_covars,
  mu,
  tau,
  calc_var = TRUE
) { 
  
  # Fit MCF.
  fit <- CalcMCF(
    time = time,
    status = status,
    idx = idx,
    calc_var = calc_var
  )
  
  # Find Area.
  area <- FindAUC(
    times = fit$time,
    values = fit$mcf,
    tau = tau
  )
  
  # Calculate influence contributions.
  # Note that psi is needed to calculate the argumentation estimator. 
  psi <- CalcPsi.AUC(
    mcf = fit,
    time = time,
    status = status,
    idx = idx,
    tau = tau
  )
  
  # Covariates.
  subj_covars_psi <- merge(
    x = subj_covars,
    y = psi,
    by = "idx"
  )
  covars <- subj_covars_psi %>% select(-c("idx", "psi")) %>% data.matrix
  
  # Find variance of area.
  var_area <- mean(subj_covars_psi$psi^2)
  
  # Resid matrix.
  n <- nrow(psi)
  center <- matrix(data = mu, nrow = n, ncol = length(mu), byrow = TRUE)
  resid <- covars - center
  
  # Covariance.
  sigma <- t(resid) %*% resid / (n^2)
  gamma <- t(resid) %*% subj_covars_psi$psi / (n^2)
  mean_resid <- t(resid) %*% rep(1, n) / n
  
  # Output.
  out <- list(
    "n" = n,
    "mcf" = fit,
    "area" = area,
    "var_area" = var_area,
    "sigma" = sigma,
    "gamma" = gamma,
    "mean_resid" = mean_resid
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Calculate Test Statistics for Augmentation Estimator.
#' 
#' Calculate test statistics for augmentation estimator.
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to remain
#'   at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param covars Per-subject covariate matrix.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param return_areas Return the AUCs?
#' @importFrom dplyr "%>%" group_by select summarise summarise_all
#' @importFrom MASS ginv
#' @export 
#' @return If `return_areas`, list containing:
#' \itemize{
#'   \item 'marg_area', marginal area for each arm.
#'   \item 'mcf', mean cumulative function for each strata.
#'   \item 'contrasts', including the augmented difference of areas.
#' }
#'  Else, only 'contrasts' is returned. 

AUC.Stats.Aug <- function(
  time,
  status,
  arm,
  idx,
  covars, 
  tau, 
  alpha, 
  return_areas = FALSE
) {
  
  # Form data.
  data <- data.frame(time, status, arm, idx)
  covars <- data.frame(arm, idx, covars)
  
  # Summarize to per-subject covariate data.
  subj_covars <- covars %>% 
    dplyr::group_by(idx) %>%
    dplyr::summarise_all(.funs = mean, .groups = "drop") %>%
    as.data.frame
  
  mu <- subj_covars %>%
    dplyr::select(-c("idx", "arm")) %>%
    dplyr::summarise_all(.funs = mean) %>%
    as.numeric
  
  # Split data.
  data1 <- subset(x = data, arm == 1)
  subj_covars1 <- subset(x = subj_covars, arm == 1, select = -arm)
  
  data0 <- subset(x = data, arm == 0)
  subj_covars0 <- subset(x = subj_covars, arm == 0, select = -arm)
  
  # Areas. 
  a1 <- AUC.Area.Aug(
    time = data1$time,
    status = data1$status,
    idx = data1$idx,
    subj_covars = subj_covars1,
    mu = mu,
    tau = tau,
    calc_var = return_areas
  )
  a0 <- AUC.Area.Aug(
    time = data0$time,
    status = data0$status,
    idx = data0$idx,
    subj_covars = subj_covars0,
    mu = mu,
    tau = tau,
    calc_var = return_areas
  )
  
  # Augmentation difference.
  sigma <- a1$sigma + a0$sigma
  gamma <- a1$gamma + a0$gamma
  omega <- MASS::ginv(sigma) %*% gamma
  correction <- as.numeric(
    t(omega) %*% (a1$mean_resid - a0$mean_resid)
  )
  delta <- a1$area - a0$area - correction
  var_delta <- a1$var_area / a1$n + 
    a0$var_area / a0$n - 
    as.numeric(t(omega) %*% gamma)
  se_delta <- sqrt(var_delta)
  
  # Contrast data.frame.
  crit <- qnorm(p = 1 - alpha / 2)
  contrasts <- data.frame(
    "contrast" = "A1-A0",
    "observed" = delta,
    "se" = se_delta,
    "lower" = delta - crit * se_delta,
    "upper" = delta + crit * se_delta,
    "p" = 2 * pnorm(q = abs(delta) / se_delta, lower.tail = FALSE)
  )
  
  # Output
  if(return_areas){
    
    # Areas.
    marg_areas <- data.frame(
      "arm" = c(0, 1),
      "n" = c(a0$n, a1$n),
      "tau" = tau,
      "area" = c(a0$area, a1$area),
      "se" = c(sqrt(a0$var_area / a0$n), sqrt(a1$var_area / a1$n))
    )
    
    # MCF.
    mcf1 <- a1$mcf
    mcf1$arm <- 1
    mcf0 <- a0$mcf
    mcf0$arm <- 0
    
    # Outputs.
    out <- list(
      "marg_areas" = marg_areas,
      "mcf" = rbind(mcf1, mcf0),
      "contrasts" = contrasts
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
#'   Note that subjects who are neither censored nor die are assumed to remain
#'   at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param tau Truncation time.
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
#' @param obs_stats Observed contrasts.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return Data.frame containing:
#' \itemize{
#'   \item Bootstrap difference 'boot_diff' and ratio 'boot_ratio' of areas.
#'   \item An indicator that the bootstrap difference was of the opposite
#'     sign, '1side_boot_diff'.
#' }

Boot.Sim.Aug <- function(
  time,
  status,
  arm,
  idx,
  covars,
  obs_stats,
  tau,
  alpha,
  reps
) {
  
  # Partition data.
  data <- data.frame(
    time,
    status,
    arm,
    idx,
    covars
  )
  data1 <- subset(x = data, arm == 1)
  data0 <- subset(x = data, arm == 0)
  n0 <- nrow(data0)
  
  # Bootstrap function.
  loop <- function(b) {
    
    # Bootstrap data sets.
    boot0 <- GroupBoot(data0)
    boot1 <- GroupBoot(data1)
    boot <- rbind(boot1, boot0)
    
    # Bootstrap statistics.
    boot_stats <- AUC.Stats.Aug(
      time = boot$time,
      status = boot$status,
      arm = boot$arm,
      idx = boot$idx,
      covars = boot %>% select(-c("time", "status", "arm", "idx")),
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
  
  sim <- lapply(seq_len(reps), loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c("boot_diff", "is_diff_sign")
  return(sim)
}


#' Bootstrap Confidence Intervals for Augmentation Estimator.
#'
#' @param sim Bootstrap simulation results, as generated by \code{\link{Boot.Sim.Aug}}.
#' @param obs_stats Observed contrasts.
#' @param alpha Type I error.
#' @importFrom stats sd
#' @return Data.framne containing the equi-talied and highest-density bootstrap
#'   confidence intervals.

Boot.CIs.Aug <- function(
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
  
  # HDI for difference.
  reps <- nrow(sim)
  hdi_diff <- HighDensCI(
    x = sim$boot_diff,
    alpha = alpha,
    min_tail_prob = 1 / reps,
    intervals = 1e3
  )
  
  # Format confidence intervals.
  cis <- data.frame(
    rbind(
      eti_diff,
      hdi_diff
    )
  )
  rownames(cis) <- NULL
  
  # Output.
  out <- data.frame(
    method = "bootstrap",
    type = c("equitailed", "highest-density"),
    contrast = rep(c("A1-A0"), times = 2),
    observed = obs_stats$observed,
    se = sd(sim$boot_diff),
    cis
  )
  return(out)
}

# -----------------------------------------------------------------------------

#' Permutation Inference
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to remain
#'   at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
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

CompAUCs.Perm.Aug <- function(
  time,
  status,
  arm,
  idx,
  covars,
  obs_stats,
  tau,
  alpha,
  reps
) {
  
  # Format data.
  data <- data.frame(
    time = time,
    status = status,
    arm = arm,
    idx = idx,
    covars = covars
  )
  
  # Permutation function.
  loop <- function(b) {
    
    # Permute data.
    perm <- PermData(data)
    
    # Permutation statistics
    perm_stats <- AUC.Stats.Aug(
      time = perm$time,
      status = perm$status,
      arm = perm$arm,
      idx = perm$idx,
      covars = perm %>% select(-c("time", "status", "arm", "idx")),
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
  
  sim <- lapply(seq_len(reps), loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c(
    "perm_diff",
    "perm_1sided"
  )
  return(sim)
}


# -----------------------------------------------------------------------------

#' Inference on the Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0) with adjusted for covariates.
#'
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death.
#'   Note that subjects who are neither censored nor die are assumed to remain
#'   at risk throughout the observation period.
#' @param arm Arm, coded as 1 for treatment, 0 for reference.
#' @param idx Unique subject index.
#' @param tau Truncation time.
#' @param covars Optional covariate matrix. Rows should correspond with the
#'   subject index `idx`. Factor and interaction terms should be expanded.
#' @param alpha Alpha level.
#' @param boot Logical, construct bootstrap confidence intervals?
#' @param perm Logical, perform permutation test?
#' @param reps Replicates for bootstrap/permutation inference.
#' @importFrom dplyr "%>%" select
#' @importFrom stats quantile
#' @importFrom methods new
#' @return Object of class CompareAugAUCs with these slots:
#' \itemize{
#'   \item `@Areas`: Marginal AUC for each arm.
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals.
#'   \item `@MCF`: Mean cumulative count curve for each arm.
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics.
#' }
#' @examples 
#' \donttest{
#' # Simulate data set.
#' data <- GenData()
#'
#' aucs <- CompareAUCs(
#'   time = data$time,
#'   status = data$status,
#'   arm = data$arm,
#'   idx = data$idx,
#'   covars = data$covar,
#'   tau = 2,
#'   boot = TRUE,
#'   perm = TRUE,
#'   reps = 25,
#'   alpha = 0.05
#' )
#' show(aucs)
#' }

CompareAugAUCs <- function(
  time,
  status,
  arm,
  idx,
  tau,
  covars = NULL,
  alpha = 0.05,
  boot = FALSE,
  perm = FALSE,
  reps = 2000
) {
  
  if (is.null(covars)) {
    covars <- rep(1, length(idx))
  }
  
  # Observed test stats.
  obs <- AUC.Stats.Aug(
    time = time,
    status = status,
    arm = arm,
    idx = idx,
    covars = covars,
    tau = tau,
    alpha = alpha,
    return_areas = TRUE
  )
  obs_stats <- obs$contrasts
  
  # CIs.
  cis <- obs_stats %>% select(-c("p"))
  cis <- data.frame(
    "method" = "asymptotic",
    "type" = "equitailed",
    cis,
    "alpha_lower" = alpha / 2,
    "alpha_upper" = alpha / 2
  )
  
  # P-values.
  pvals <- obs_stats %>% select(c("contrast", "observed", "p"))
  pvals <- cbind(
    "method" = "asymptotic",
    pvals
  )
  
  # Simulation replicates.
  sim_reps <- list()
  
  # -------------------------------------------------------
  
  # Bootstrap inference.
  if (boot) {
    
    # Simulate.
    boot_sim <- Boot.Sim.Aug(
      time = time,
      status = status,
      arm = arm,
      idx = idx,
      covars = covars,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps
    )
    sim_reps$boot_sim <- boot_sim
    
    # Confidence intervals.
    boot_cis <- Boot.CIs.Aug(
      sim = boot_sim,
      obs_stats = obs_stats,
      alpha = alpha
    )
    cis <- rbind(
      cis,
      boot_cis
    )
    cis <- cis[order(cis$contrast), ]
    
    # P-value.
    boot_p <- min(2  * mean(c(1, boot_sim$is_diff_sign)), 1)
    boot_pvals <- data.frame(
      "method" = "bootstrap",
      pvals %>% select("contrast", "observed"),
      "p" = boot_p
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
    perm_sim <- CompAUCs.Perm.Aug(
      time = time,
      status = status,
      arm = arm,
      idx = idx,
      covars = covars,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps
    )
    sim_reps$perm_sim <- perm_sim
    
    # Permutation p-values.
    perm_pval <- min(2 * mean(c(1, perm_sim$perm_1sided)), 1)
    perm_pvals <- data.frame(
      "method" = "permutation",
      "contrast" = "A1-A0",
      "observed" = obs_stats$observed[1],
      "p" = perm_pval
    )
    pvals <- rbind(
      pvals,
      perm_pvals
    )
  }
  
  # -------------------------------------------------------
  
  # Output
  out <- new(
    Class = "CompAugAUCs",
    Areas = obs$marg_areas,
    CIs = cis,
    MCF = obs$mcf,
    Reps = sim_reps,
    Pvals = pvals
  )
}
