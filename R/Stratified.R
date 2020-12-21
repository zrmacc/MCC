# Purpose: Stratified estimation of AUC.
# Updated: 2020-12-17

# -----------------------------------------------------------------------------
# Calculate Stratified AUC.
# -----------------------------------------------------------------------------

#' Calculate Test Statistics for Stratified Estimator.
#'
#' @param data Data.frame containing: idx, time, status, arm, strata.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param return_areas Return the AUCs?
#' @importFrom dplyr "%>%" filter group_by group_split inner_join select
#'   summarise
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

CalcStratAUC <- function(
  data,
  tau,
  alpha,
  return_areas = FALSE
) {
  
  # Stratum sizes.
  idx <- time <- status <- arm <- strata <- NULL
  stratum_sizes <- data %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise("n" = length(unique(idx)), .groups = "drop") %>%
    as.data.frame
  stratum_sizes$weight <- stratum_sizes$n / sum(stratum_sizes$n)
  
  # Stratum areas.
  stratum_areas <- data %>%
    dplyr::group_by(arm, strata) %>%
    dplyr::summarise(
      StratumAUC(time, status, idx, tau = tau, calc_var = return_areas),
      .groups = "drop"
    ) %>% 
    dplyr::inner_join(
      stratum_sizes[, c("strata", "weight")], 
      by = "strata"
    ) %>%
    as.data.frame
  
  # Marginal areas.
  area <- n <- se_area <- weight <- NULL
  marg_areas <- stratum_areas %>%
    dplyr::group_by(arm) %>%
    dplyr::summarise(
      MargAUC(areas = area, ses = se_area, weights = weight, n = n),
      .groups = "drop"
    ) %>%
    as.data.frame
  marg_areas$tau <- tau
  
  # Difference and ratio.
  contrasts <- ContrastAreas(
    marg_areas = marg_areas,
    alpha = alpha
  )
  
  # Output
  if(return_areas){
    
    # Marginal MCF for arm 1.
    mcf1 <- data %>%
      dplyr::filter(arm == 1) %>%
      dplyr::group_by(strata) %>%
      dplyr::summarise(
        CalcMCF(time, status, idx),
        .groups = "keep"
      ) %>%
      dplyr::group_split()
    avg_mcf1 <- AvgMCF(mcf1, weights = stratum_sizes$weight)
    avg_mcf1$arm <- 1
    
    # Marginal MCF for arm 0.
    mcf0 <- data %>%
      dplyr::filter(arm == 0) %>%
      dplyr::group_by(strata) %>%
      dplyr::summarise(
        CalcMCF(time, status, idx),
        .groups = "keep"
      ) %>%
      dplyr::group_split()
    avg_mcf0 <- AvgMCF(mcf1, weights = stratum_sizes$weight)
    avg_mcf0$arm <- 0
    
    avg_mcf <- rbind(avg_mcf1, avg_mcf0)
    
    # Outputs.
    out <- list(
      avg_mcf = avg_mcf,
      contrasts = contrasts,
      marg_areas = marg_areas,
      stratum_areas = stratum_areas
    )
  } else {
    out <- contrasts
  }
  return(out)
}


# -----------------------------------------------------------------------------

#' Stratum AUC.
#' 
#' Calculates the AUC for a single stratum.
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#' @param idx Unique subject index. 
#' @param tau Truncation time. 
#' @param calc_var Calculate analytical variance? 
#' @return Data.frame containing:
#' \itemize{
#'   \item Truncation time 'tau' and estimated 'area'.
#'   \item Variance 'var_area' and standard error 'se_area', if requested.
#' }

StratumAUC <- function(
  time, 
  status,
  idx,
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
  
  # Find AUC.
  area <- AUC(
    times = fit$time,
    values = fit$mcf,
    tau = tau
  )
  
  # Output.
  n <- length(unique(idx))
  out <- data.frame(
    "n" = n,
    "tau" = tau,
    "area" = area
  )
  
  if (calc_var) {
    
    # Calculate influence contributions.
    psi <- data.frame(
      time = time,
      status = status,
      idx = idx
    ) %>%
      dplyr::group_by(idx) %>%
      dplyr::summarise(
        "psi" = PsiAUCi(mcf = fit, time = time, status = status, tau = tau),
        .groups = "drop" ) %>%
      dplyr::pull(psi) 
    
    # Find variance of area.
    out$var_area <- mean(psi^2)
    out$se_area <- sqrt(out$var_area / n)
  }
  
  return(out)
}


# -----------------------------------------------------------------------------

#' Calculate Marginal Area
#' 
#' Calculates the marginal AUC across strata.
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
  
  # Output.
  out <- data.frame(
    "n" = sum(n),
    "area" = sum(weights * areas)
  )
  
  # Add standard error, if provided.
  if (!is.null(ses)) {
    out$se <- sqrt(sum(weights^2 * ses^2))
  }
  
  return(out)
}


# -----------------------------------------------------------------------------
# Main function for stratified estimator.
# -----------------------------------------------------------------------------

#' Inference on the Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0).
#'
#' @param data Formatted data.frame. \code{\link{FormatData}}.
#' @param tau Truncation time.
#' @param alpha Alpha level.
#' @param boot Logical, construct bootstrap confidence intervals?
#' @param perm Logical, perform permutation test?
#' @param reps Replicates for bootstrap/permutation inference.
#' @importFrom stats quantile
#' @importFrom methods new
#' @importFrom dplyr "%>%" select
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
#' aucs <- CompareAUCs(
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
  data,
  tau,
  alpha = 0.05,
  boot = FALSE,
  perm = FALSE,
  reps = 2000
) {

  # Observed test stats.
  obs <- CalcStratAUC(
    data = data,
    alpha = alpha,
    tau = tau,
    return_areas = TRUE
  )
  obs_stats <- obs$contrasts

  # CIs.
  cis <- obs_stats %>% dplyr::select(-c("p"))
  cis <- data.frame(
    method = "asymptotic",
    cis
  )

  # P-values.
  pvals <- obs_stats %>% dplyr::select(c("contrast", "observed", "p"))
  pvals <- data.frame(
    "method" = "asymptotic",
    pvals
  )

  # Simulation replicates.
  sim_reps <- list()

  # -------------------------------------------------------

  # Bootstrap inference.
  if (boot) {

    # Simulate.
    boot_sim <- BootSimStrat(
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
    boot_pvals <- cbind(
      "method" = "bootstrap",
      pvals %>% dplyr::select(c("contrast", "observed")),
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
    perm_sim <- PermSimStrat(
      data = data,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps
    )
    sim_reps$perm_sim <- perm_sim

    # Permutation p-values.
    perm_pvals <- perm_sim %>%
      dplyr::select("perm_diff_1sided", "perm_ratio_1sided") %>%
      dplyr::summarise_all(CalcP) %>% as.numeric
      
    perm_pvals <- data.frame(
      "method" = "permutation",
      "contrast" = c("A1-A0", "A1/A0"),
      "observed" = obs_stats$observed,
      "p" = perm_pvals
    )
    rownames(perm_pvals) <- NULL
    pvals <- rbind(
      pvals,
      perm_pvals
    )
    pvals <- pvals[order(pvals$observed), ]
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
