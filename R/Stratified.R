# Purpose: Stratified estimation of AUC.
# Updated: 2022-05-19

# -----------------------------------------------------------------------------
# Calculate Stratified AUC.
# -----------------------------------------------------------------------------

#' Calculate Test Statistics for Stratified Estimator
#'
#' @param data Data.frame containing: {arm, idx, status, strata, time}.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param return_areas Return the AUCs?
#' @return If `return_areas`, list containing:
#' \itemize{
#'   \item 'avg_mcf', average MCF across strata.
#'   \item 'contrasts', including the difference and ratio of areas.
#'   \item 'marg_areas', the marginal areas for each arm.
#'   \item 'stratum_areas', the per-stratum and per-arm areas.
#'   \item 'weights', the weights of each stratum.
#' }
#' @importFrom dplyr "%>%"
#' @noRd
CalcStratAUC <- function(
  data,
  tau,
  alpha = 0.05,
  return_areas = FALSE
) {
  
  # Stratum sizes.
  idx <- time <- status <- arm <- strata <- NULL
  stratum_sizes <- data %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise("n" = length(unique(idx)), .groups = "drop") %>%
    as.data.frame()
  stratum_sizes$weight <- stratum_sizes$n / sum(stratum_sizes$n)
  
  # Stratum areas.
  stratum_areas <- data %>%
    dplyr::group_by(arm, strata) %>%
    dplyr::summarise(
      StratumAUC(
        idx = idx,
        status = status,
        time = time,
        tau = tau,
        calc_var = return_areas
      ),
      .groups = "drop"
    ) %>% 
    dplyr::inner_join(
      stratum_sizes[, c("strata", "weight")], 
      by = "strata"
    ) %>%
    as.data.frame()
  
  # Marginal areas.
  area <- n <- se_area <- weight <- NULL
  marg_areas <- stratum_areas %>%
    dplyr::group_by(arm) %>%
    dplyr::summarise(
      MargAUC(areas = area, n = n, ses = se_area, weights = weight),
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
  if (return_areas) {
    
    avg_mcf <- CalcMargMCF(data)
    
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
#' @param idx Subject index.
#' @param status Event status.
#' @param time Observation time.
#' @param tau Truncation time. 
#' @param calc_var Calculate analytical variance? 
#' @return Data.frame containing:
#' \itemize{
#'   \item Truncation time 'tau' and estimated 'area'.
#'   \item Variance 'var_area' and standard error 'se_area', if requested.
#' }
#' @noRd

StratumAUC <- function(
  idx,
  status,
  time,
  tau,
  calc_var = TRUE
) { 
  
  # Fit MCF.
  mcf <- CalcMCF(
    idx = idx,
    status = status,
    time = time,
    calc_var = calc_var
  )
  
  # Find AUC.
  area <- AUC(
    times = mcf$time,
    values = mcf$mcf,
    tau = tau
  )
  
  # Output.
  n <- length(unique(idx))
  out <- data.frame(
    n = n,
    tau = tau,
    area = area
  )
  
  if (calc_var) {
    
    # Find variance of area.
    df <- data.frame(idx = idx, status = status, time = time)
    out$var_area <- VarAUC(df, tau, mcf = mcf)
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
#' @param n Sample sizes.
#' @param ses Standard errors.
#' @param weights Weights.
#' @return Data.frame containing:
#' \itemize{
#'   \item The marginal area.
#'   \item The standard error of the area.
#' }
#' @noRd

MargAUC <- function(
  areas,
  n,
  ses,
  weights
) {
  
  # Output.
  out <- data.frame(n = sum(n), area = sum(weights * areas))
  
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
#' @importFrom dplyr "%>%" 
#' @examples
#' \donttest{
#' # Simulate data set.
#' n <- 100
#' covariates <- data.frame(
#'   arm = c(rep(1, n/2), rep(0, n/2)),
#'   strata = stats::rbinom(n, 1, 0.2)
#' )
#' data <- GenData(
#'   beta_event = c(log(0.5), log(0.8)),
#'   covariates = covariates
#' )
#'
#' aucs <- CompareAUCs(
#'   data,
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
  cis <- data.frame(method = "asymptotic", cis)

  # P-values.
  pvals <- obs_stats %>% dplyr::select(c("contrast", "observed", "p"))
  pvals <- data.frame(method = "asymptotic", pvals)

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
      method = "bootstrap",
      pvals %>% dplyr::select(c("contrast", "observed")),
      p = boot_p
    )
    pvals <- rbind(pvals, boot_pvals)
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
      method = "permutation",
      contrast = c("A1-A0", "A1/A0"),
      observed = obs_stats$observed,
      p = perm_pvals
    )
    rownames(perm_pvals) <- NULL
    pvals <- rbind(pvals, perm_pvals)
    pvals <- pvals[order(pvals$observed), ]
  }

  # -------------------------------------------------------

  # Output
  out <- methods::new(
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
