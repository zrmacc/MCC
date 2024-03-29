# Purpose: Augmentation estimator of AUC.
# Updated: 2024-02-19

# -----------------------------------------------------------------------------
# Calculated Augmented AUC.
# -----------------------------------------------------------------------------

#' Calculate Test Statistics for Augmentation Estimator.
#' 
#' Calculate test statistics for augmentation estimator.
#'
#' @param data Data.frame containing {arm, idx, status, time, weights}.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param return_areas Return the AUCs?
#' @return If return_areas is TRUE, list containing:
#' \itemize{
#'   \item 'marg_area', marginal area for each arm.
#'   \item 'mcf', mean cumulative function for each strata.
#'   \item 'contrasts', including the augmented difference of areas.
#' }
#' @export 
CalcAugAUC <- function(
  data,
  tau, 
  alpha = 0.05, 
  return_areas = FALSE
) {
  
  # Summarize to per-subject covariate data.
  arm <- idx <- NULL
  covars <- data %>% 
    dplyr::select(-c("time", "status", "weights")) %>%
    dplyr::group_by(arm, idx) %>%
    dplyr::summarise_all(.funs = mean, .groups = "drop") %>%
    as.data.frame()
  
  mu <- covars %>%
    dplyr::select(-c("arm", "idx")) %>%
    dplyr::summarise_all(.funs = mean) %>%
    as.numeric()
  
  # Areas. 
  a1 <- AugAUC(
    data = data %>% dplyr::filter(arm == 1),
    covars = covars %>% dplyr::filter(arm == 1) %>% dplyr::select(-arm),
    mu = mu,
    tau = tau,
    calc_var = return_areas
  )
  a0 <- AugAUC(
    data = data %>% dplyr::filter(arm == 0),
    covars = covars %>% dplyr::filter(arm == 0) %>% dplyr::select(-arm),
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
  crit <- stats::qnorm(p = 1 - alpha / 2)
  contrasts <- data.frame(
    contrast = "A1-A0",
    observed = delta,
    se = se_delta,
    lower = delta - crit * se_delta,
    upper = delta + crit * se_delta,
    p = 2 * stats::pnorm(q = abs(delta) / se_delta, lower.tail = FALSE)
  )
  
  # Output
  if(return_areas){
    
    # Areas.
    marg_areas <- data.frame(
      arm = c(0, 1),
      n = c(a0$n, a1$n),
      tau = tau,
      area = c(a0$area, a1$area),
      se = c(sqrt(a0$var_area / a0$n), sqrt(a1$var_area / a1$n))
    )
    
    # MCF.
    mcf1 <- a1$mcf
    mcf1$arm <- 1
    mcf0 <- a0$mcf
    mcf0$arm <- 0
    
    # Outputs.
    out <- list(
      marg_areas = marg_areas,
      mcf = rbind(mcf1, mcf0),
      contrasts = contrasts
    )
  } else {
    out <- contrasts
  }
  return(out)
}


# -----------------------------------------------------------------------------

#' Calculate Augmentation Components for a Single Arm
#'
#' @param data Data.frame including {idx, status, time, weights}.
#' @param covars Per-subject covariate data.
#' @param mu Grand mean.
#' @param tau Truncation time. 
#' @param calc_var Calculate analytical variance of MCF?
#' @return List containing:
#' \itemize{
#'   \item Tabulated mean cumulative function 'mcf'.
#'   \item Estimated 'area' and its variance 'var_area'.
#'   \item 'sigma', \eqn{\frac{1}{n^2}\sum_{i=1}^{n}(X_{i}-\mu)(X_{i}-\mu)'}.
#'   \item 'gamma', \eqn{\frac{1}{n^2}\sum_{i=1}^{n}(X_{i}-\mu)\xi_{i}}.
#'   \item 'mean_resid', \eqn{\frac{1}{n}\sum_{i=1}^{n}(X_{i}-\mu)}.
#' }
#' @noRd
AugAUC <- function(
  data,
  covars,
  mu,
  tau,
  calc_var = TRUE
) { 
  
  # Fit MCF.
  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    weights = data$weights,
    calc_var = calc_var
  )
  
  # Find Area.
  area <- AUC(
    times = mcf$time,
    values = mcf$mcf,
    tau = tau
  )
  
  # Calculate influence contributions.
  psi <- VarAUC(data, tau, mcf = mcf, return_psi = TRUE)
  
  # Covariates.
  covars <- covars %>% 
    dplyr::inner_join(psi, by = "idx") %>%
    dplyr::select(-c("idx", "psi")) %>%
    data.matrix()
  
  # Find variance of area.
  var_area <- mean(psi$psi^2)
  
  # Calculate remaining augmentation components.
  aug_comp <- CalcAugComp(covars, mu, psi$psi)
  
  # Output.
  out <- list(
    area = area,
    gamma = aug_comp$gamma,
    n = aug_comp$n,
    mean_resid = aug_comp$mean_resid,
    mcf = mcf,
    sigma = aug_comp$sigma,
    var_area = var_area
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
#' @param data Formatted data frame. \code{\link{FormatData}}.
#' @param tau Truncation time.
#' @param alpha Alpha level.
#' @param boot Logical, construct bootstrap confidence intervals?
#' @param perm Logical, perform permutation test?
#' @param reps Replicates for bootstrap/permutation inference.
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
#' n <- 100
#' covariates <- data.frame(
#'   arm = c(rep(1, n/2), rep(0, n/2)),
#'   covar = rnorm(n)
#' )
#' data <- GenData(
#'   beta_event = c(log(0.5), 1),
#'   covariates = covariates
#' )
#'
#' aucs <- CompareAUCs(
#'   data,
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
  data,
  tau,
  alpha = 0.05,
  boot = FALSE,
  perm = FALSE,
  reps = 2000
) {
  
  # Observed test stats.
  obs <- CalcAugAUC(
    data = data,
    tau = tau,
    alpha = alpha,
    return_areas = TRUE
  )
  obs_stats <- obs$contrasts
  
  # CIs.
  cis <- obs_stats %>% dplyr::select(-c("p"))
  cis <- data.frame(method = "asymptotic", cis)
  
  # P-values.
  pvals <- obs_stats %>% dplyr::select(c("contrast", "observed", "p"))
  pvals <- cbind(method = "asymptotic", pvals)
  
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
      method = "bootstrap",
      pvals %>% dplyr::select("contrast", "observed"),
      p = boot_p
    )
    pvals <- rbind(pvals, boot_pvals)
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
      method = "permutation",
      contrast = "A1-A0",
      observed = obs_stats$observed[1],
      p = perm_pval
    )
    pvals <- rbind(
      pvals,
      perm_pvals
    )
  }
  
  # -------------------------------------------------------
  
  # Output
  out <- methods::new(
    Class = "CompAugAUCs",
    Areas = obs$marg_areas,
    CIs = cis,
    MCF = obs$mcf,
    Reps = sim_reps,
    Pvals = pvals
  )
}
