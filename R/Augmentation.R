# Purpose: Augmentation estimator of AUC.
# Updated: 2020-12-17

# -----------------------------------------------------------------------------
# Calculated Augmented AUC.
# -----------------------------------------------------------------------------

#' Calculate Test Statistics for Augmentation Estimator.
#' 
#' Calculate test statistics for augmentation estimator.
#'
#' @param data Data.frame containing: idx, time, status, arm, and covars.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param return_areas Return the AUCs?
#' @importFrom dplyr "%>%" group_by select summarise summarise_all
#' @importFrom MASS ginv
#' @export 
#' @return If return_areas is TRUE, list containing:
#' \itemize{
#'   \item 'marg_area', marginal area for each arm.
#'   \item 'mcf', mean cumulative function for each strata.
#'   \item 'contrasts', including the augmented difference of areas.
#' }
#'  Else, only 'contrasts' is returned. 

CalcAugAUC <- function(
  data,
  tau, 
  alpha, 
  return_areas = FALSE
) {
  
  # Summarize to per-subject covariate data.
  arm <- idx <- NULL
  subj_covars <- data %>% 
    dplyr::select(-c("time", "status")) %>%
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
  a1 <- AugAUC(
    data = data1,
    subj_covars = subj_covars1,
    mu = mu,
    tau = tau,
    calc_var = return_areas
  )
  a0 <- AugAUC(
    data = data0,
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

#' AUC Calculation for a Single Arm.
#'
#' @param data Data.frame including: time, status, idx.
#' @param subj_covars Per-subject covariate data.
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

AugAUC <- function(
  data,
  subj_covars,
  mu,
  tau,
  calc_var = TRUE
) { 
  
  # Fit MCF.
  fit <- CalcMCF(
    time = data$time,
    status = data$status,
    idx = data$idx,
    calc_var = calc_var
  )
  
  # Find Area.
  area <- AUC(
    times = fit$time,
    values = fit$mcf,
    tau = tau
  )
  
  # Calculate influence contributions.
  # Note that psi is needed to calculate the argumentation estimator. 
  idx <- time <- status <- NULL
  subj_covars_psi <- data %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      "psi" = PsiAUCi(mcf = fit, time = time, status = status, tau = tau),
      .groups = "drop" ) %>%
    dplyr::inner_join(y = subj_covars, by = "idx")
  
  # Covariates.
  covars <- subj_covars_psi %>% dplyr::select(-c("idx", "psi")) %>% data.matrix
  
  # Find variance of area.
  var_area <- mean(subj_covars_psi$psi^2)
  
  # Resid matrix.
  n <- nrow(subj_covars)
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
# Main function for aumentation estimator.
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
  
  # Create data.frames.
  data <- data.frame(
    idx = idx,
    time = time,
    status = status,
    arm = arm,
    covars
  )
  
  # Observed test stats.
  obs <- CalcAugAUC(
    data = data,
    tau = tau,
    alpha = alpha,
    return_areas = TRUE
  )
  obs_stats <- obs$contrasts
  
  # CIs.
  cis <- obs_stats %>% select(-c("p"))
  cis <- data.frame(
    "method" = "asymptotic",
    cis
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
    boot_sim <- BootSimAug(
      data = data,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps
    )
    sim_reps$boot_sim <- boot_sim
    
    # Confidence intervals.
    boot_cis <- BootCIs(
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
    boot_p <- CalcP(boot_sim$is_diff_sign)
    boot_pvals <- data.frame(
      "method" = "bootstrap",
      pvals %>% dplyr::select("contrast", "observed"),
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
    perm_sim <- PermSimAug(
      data = data,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps
    )
    sim_reps$perm_sim <- perm_sim
    
    # Permutation p-values.
    perm_pval <- CalcP(perm_sim$perm_1sided)
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
