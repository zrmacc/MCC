# Purpose: Stratified estimation of AUC.
# Updated: 2024-02-20

# -----------------------------------------------------------------------------
# Calculate Stratified AUC.
# -----------------------------------------------------------------------------

#' Calculate Test Statistics for Stratified Estimator
#'
#' @param data Data.frame containing: {arm, idx, status, strata, time, weights}.
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
#' @noRd
CalcStratAUC <- function(
  data,
  tau,
  alpha = 0.05,
  return_areas = FALSE
) {
  
  # Stratum sizes.
  idx <- time <- status <- arm <- strata <- weights <- NULL
  stratum_sizes <- data %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise("n" = length(unique(idx)), .groups = "drop") %>%
    as.data.frame()
  stratum_sizes$strat_weight <- stratum_sizes$n / sum(stratum_sizes$n)
  
  # Stratum areas.
  stratum_areas <- data %>%
    dplyr::group_by(arm, strata) %>%
    dplyr::reframe(
      StratumAUC(idx, status, time, tau, weights)
    ) %>% 
    dplyr::inner_join(
      stratum_sizes[, c("strata", "strat_weight")], 
      by = "strata"
    ) %>%
    as.data.frame()
  
  # Marginal areas.
  area <- n <- se_area <- strat_weight <- NULL
  marg_areas <- stratum_areas %>%
    dplyr::group_by(arm) %>%
    dplyr::reframe(
      MargAUC(areas = area, n = n, ses = se_area, weights = strat_weight),
    ) %>%
    as.data.frame()
  marg_areas$tau <- tau
  
  # Difference and ratio.
  contrasts <- ContrastAreas(marg_areas = marg_areas, alpha = alpha) 
  
  # Output
  if (return_areas) {
    
    if (nrow(marg_areas) == 2) {
      avg_mcf <- CalcMargMCF(data)
    } else {
      avg_mcf <- CalcMCF(
        idx = data$idx,
        status = data$status,
        time = data$time,
        weights = data$weights
      )
    }
    
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
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
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
    weights,
    calc_var = TRUE
) { 
  
  # Fit MCF.
  mcf <- CalcMCF(
    idx = idx,
    status = status,
    time = time,
    weights = weights,
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
    out$var_area <- VarAUC(df, tau, mcf = mcf, weights = weights)
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
#' @param weights Stratum size weights.
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
# Average MCF and marginal MCF (from AvgMCF.R)
# -----------------------------------------------------------------------------

#' Calculate Average MCF Curve
#'
#' Calculates the weighted average of MCF curves for a stratified analysis.
#'
#' @param curve_list List of tabulated MCFs as returned by \code{\link{CalcMCF}}.
#' @param strat_weights Numeric vector of stratum weights.
#' @return Data.frame containing `Time` and the averaged MCF `Avg_MCF`.
#' @noRd
AvgMCF <- function(curve_list, strat_weights) {
  time <- lapply(curve_list, function(x) { x$time })
  time <- do.call(c, time)
  time <- sort(unique(time))
  aux <- function(x) {
    g <- stats::stepfun(x$time, c(0, x$mcf), right = FALSE)
    return(g(time))
  }
  mcfs <- lapply(curve_list, aux)
  mcfs <- do.call(cbind, mcfs)
  avg_mcf <- mcfs %*% strat_weights
  aux <- function(x) {
    g <- stats::stepfun(x$time, c(0, x$var_mcf), right = FALSE)
    return(g(time))
  }
  vars <- lapply(curve_list, aux)
  vars <- do.call(cbind, vars)
  avg_var <- vars %*% (strat_weights^2)
  out <- data.frame(
    time = time,
    mcf = avg_mcf,
    var_mcf = avg_var
  )
  out$se_mcf <- sqrt(out$var_mcf)
  return(out)
}

#' Calculate Marginal MCF
#'
#' Calculates the marginal MCF, averaged across strata, with stratum
#' weights proportional to the total number of subjects (across arms)
#' belonging to that stratum.
#'
#' @param data Data.frame containing (arm, idx, status, strata, time, weights).
#' @return Data.frame.
#' @export
CalcMargMCF <- function(data) {
  arm <- idx <- status <- strata <- time <- weights <- NULL
  data <- data %>%
    dplyr::select(arm, idx, status, strata, time, weights)
  n0 <- n1 <- n <- w <- NULL
  stratum_sizes <- data %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      n0 = length(unique(idx[arm == 0])),
      n1 = length(unique(idx[arm == 1])),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n = n0 + n1,
      w = n / sum(n),
      w0 = w / sum(w[n0 != 0]),
      w1 = w / sum(w[n1 != 0])
    )
  stratum_sizes$w1[stratum_sizes$n1 == 0] <- 0
  stratum_sizes$w0[stratum_sizes$n0 == 0] <- 0
  mcf1 <- data %>%
    dplyr::filter(arm == 1) %>%
    dplyr::group_by(strata) %>%
    dplyr::reframe(
      CalcMCF(idx = idx, status = status, time = time, weights = weights, calc_var = TRUE)
    ) %>%
    dplyr::group_by(strata) %>%
    dplyr::group_split()
  avg_mcf1 <- AvgMCF(mcf1, strat_weights = stratum_sizes$w1[stratum_sizes$w1 != 0])
  avg_mcf1$arm <- 1
  mcf0 <- data %>%
    dplyr::filter(arm == 0) %>%
    dplyr::group_by(strata) %>%
    dplyr::reframe(
      CalcMCF(idx = idx, status = status, time = time, weights = weights, calc_var = TRUE)
    ) %>%
    dplyr::group_by(strata) %>%
    dplyr::group_split()
  avg_mcf0 <- AvgMCF(mcf0, strat_weights = stratum_sizes$w0[stratum_sizes$w0 != 0])
  avg_mcf0$arm <- 0
  avg_mcf <- rbind(avg_mcf1, avg_mcf0)
  return(avg_mcf)
}


# -----------------------------------------------------------------------------
# Main function for single-arm stratified estimator.
# -----------------------------------------------------------------------------

#' Single-Arm Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0).
#'
#' @param data Formatted data.frame. \code{\link{FormatData}}.
#' @param tau Truncation time.
#' @param alpha Alpha level.
#' @param boot Logical, construct bootstrap confidence intervals?
#' @param reps Replicates for bootstrap/permutation inference.
#' @return Object of class \code{CompStratAUCs} with these slots:
#' \itemize{
#'   \item `@Areas`: The AUC for each arm.
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals.
#'   \item `@Curves`: Mean cumulative count curve for each arm; averaged across strata
#'     if present.
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics.
#'   \item `@Weights`: Per-stratum weights and AUCs.
#' }
SingleArmStratAUC <- function(
    data,
    tau,
    alpha = 0.05,
    boot = FALSE,
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


# -----------------------------------------------------------------------------
# Main function for two-arm stratified estimator.
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
#' @return Object of class \code{CompStratAUCs} with these slots:
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
  obs <- CalcStratAUC(
    data = data,
    alpha = alpha,
    tau = tau,
    return_areas = TRUE
  )
  obs_stats <- obs$contrasts

  format_perm_pvals_strat <- function(perm_sim, obs_stats) {
    perm_pvals <- perm_sim %>%
      dplyr::select("perm_diff_1sided", "perm_ratio_1sided") %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), CalcP)) %>%
      as.numeric()
    out <- data.frame(
      method = "permutation",
      contrast = c("A1-A0", "A1/A0"),
      observed = obs_stats$observed,
      p = perm_pvals
    )
    rownames(out) <- NULL
    out
  }

  res <- .AddResamplingResults(
    obs_stats = obs_stats,
    data = data,
    tau = tau,
    alpha = alpha,
    boot = boot,
    perm = perm,
    reps = reps,
    boot_sim_fun = BootSimStrat,
    perm_sim_fun = PermSimStrat,
    format_perm_pvals = format_perm_pvals_strat
  )
  pvals <- res$pvals
  pvals <- pvals[order(pvals$observed), ]

  methods::new(
    Class = "CompStratAUCs",
    StratumAreas = obs$stratum_areas,
    MargAreas = obs$marg_areas,
    CIs = res$cis,
    MCF = obs$avg_mcf,
    Reps = res$sim_reps,
    Pvals = pvals
  )
}
