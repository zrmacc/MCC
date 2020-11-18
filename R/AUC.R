# Purpose: Functions to perform inference on the difference or ratio 
# in AUCs between two mean cumulative count curves. 

# -----------------------------------------------------------------------------
# Area under the mean cumulative function.
# -----------------------------------------------------------------------------

#' Find Area Under the Mean Cumulative Count Curve. 
#'
#' @param times MCF times.
#' @param values MCF values
#' @param tau Truncation time.
#' @importFrom stats integrate stepfun 
#' @return Numeric area under the curve. 

FindAUC <- function(times, values, tau) {
  g <- stepfun(
    x = times,
    y = c(0, values)
  )
  area <- integrate(f = g, lower = 0, upper = tau, subdivisions = 2e3)
  return(area$value)
}

# -----------------------------------------------------------------------------

#' Calculate AUC
#' 
#' Used to find the AU MCF for each stratum.
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period. 
#' @param idx Unique subject index. 
#' @param tau Truncation time. 
#' @return Data.frame containing:
#' \itemize{
#'   \item Truncation time 'tau'.
#'   \item Estimated 'area', its variance 'var_area' and standard error 'se_area'.
#' }

AUC.Area <- function(
  time, 
  status,
  idx,
  tau
) { 
  
  # Fit MCF.
  fit <- CalcMCF(
    time = time,
    status = status,
    idx = idx
  )
  
  # Find Area.
  area <- FindAUC(
    times = fit$time,
    values = fit$mcf,
    tau = tau
  )
  
  # Find variance of area.
  var_area <- CalcVarAUC(
    mcf = fit,
    time = time,
    status = status,
    idx = idx,
    tau = tau
  )
  
  # Output.
  out <- data.frame(
    'tau' = tau,
    'area' = area,
    'var_area' = var_area
  )
  out$se_area <- sqrt(out$var_area / length(unique(idx)))
  return(out)
}

# -----------------------------------------------------------------------------
# Variance calculation for area under the mean cumulative function.
# -----------------------------------------------------------------------------

#' Calculate Variance Contribution of an Individual Subject to Area Under
#'   the Mean Cumulative Function.
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times for a single subject.
#' @param status Status indicator for a single subject.
#' @param tau Truncation time.
#' @return Numeric square influence function for a single observation.

CalcVarAUC.i <- function(
  mcf,
  time,
  status,
  tau
) { 

  # Subject-specific martingales.
  out <- CalcMartigales.i(
    mcf = mcf,
    time = time,
    status = status
  )
  
  # Truncate.
  out <- out[out$time <= tau, ]
  mcf <- mcf[mcf$time <= tau, ]
  
  # Calculate variance.
  int_1 <- sum((tau - out$time) * mcf$surv / mcf$prop_at_risk * out$dM_event)
  int_2 <- sum((tau - out$time) * mcf$mcf  / mcf$prop_at_risk * out$dM_death)
  psi <- int_1 + int_2
  psi2 <- psi^2
  
  # Output.
  return(psi2)
}


#' Calculate Variance of the Area Under the Mean Cumulative Function
#' 
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}.
#' @param time Observation times.
#' @param status Status indicators.
#' @param idx Unique subject index. 
#' @param tau Truncation time.
#' @return Numeric estimated variance.

CalcVarAUC <- function(
  mcf,
  time,
  status,
  idx,
  tau
) {
  
  # Partition data by subject.
  subj <- data.frame(
    "time" = time,
    "status" = status,
    "idx" = idx
  )
  subj_split <- split(subj, subj$idx)
  
  # Calculate psi squared; Ghosh and Lin (2.2).
  psi2 <- sapply(subj_split, FUN = function(df) {
    return(CalcVarAUC.i(mcf = mcf, time = df$time, status = df$status, tau = tau))
  })
  var <- mean(psi2)
  
  # Output.
  return(var)
}

# -----------------------------------------------------------------------------
# Marginal AUC.
# -----------------------------------------------------------------------------

#' Calculate Marginal Area
#' 
#' @param areas Estimated statistics.
#' @param ses Standard errors.
#' @param weights Weights.
#' @return Data.frame containing:
#' \itemize{
#'   \item The marginal area.
#'   \item The standard error of the area.
#' }

MargAUC <- function(
  areas,
  ses,
  weights
) {
  out <- data.frame(
    "area" = sum(weights * areas),
    "se_area" = sqrt(sum(weights^2 * ses^2))
  )
  return(out)
}


# -----------------------------------------------------------------------------
# Contrast AUCs.
# -----------------------------------------------------------------------------

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
    "contrast" = c("A1-A0", "A1/A0"),
    "est" = c(delta, rho),
    "lower" = c(delta_lower, rho_lower),
    "upper" = c(delta_upper, rho_upper),
    "p" = c(delta_p, rho_p)
  )
  return(out)
}


# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------

#' Calculate Test Statistics
#'
#' @param data0 Data.frame containing `time`, `status`, `idx`, `strata` for arm 0.
#' @param data1 Data.frame containing `time`, `status`, `idx`, `strata` for arm 1.
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
#'  Else, only 'contrasts' is returned. 

AUC.Stats <- function(data1, data0, tau, alpha, return_areas = FALSE) {
  
  # Partition by strata.
  data1_strata <- split(data1, data1$strata, drop = TRUE)
  data0_strata <- split(data0, data0$strata, drop = TRUE)
  
  # Stratum sizes.
  aux1 <- function(x){length(unique(x$idx))}
  data1_sizes <- sapply(data1_strata, aux1)
  data0_sizes <- sapply(data0_strata, aux1)
  sizes <- data1_sizes + data0_sizes
  weights <- sizes / sum(sizes)
  
  # Stratum areas.
  aux2 <- function(x){AUC.Area(x$time, x$status, x$idx, tau)}
  data1_areas <- lapply(data1_strata, aux2)
  data1_areas <- data.frame(do.call(rbind, data1_areas))
  data1_areas <- cbind(
    "arm" = 1,
    "stratum" = names(data1_strata),
    "n" = data1_sizes,
    data1_areas
  )
  data0_areas <- lapply(data0_strata, aux2)
  data0_areas <- data.frame(do.call(rbind, data0_areas))
  data0_areas <- cbind(
    "arm" = 0,
    "stratum" = names(data0_strata),
    "n" = data0_sizes,
    data0_areas
  )
  stratum_areas <- rbind(data1_areas, data0_areas)
  rownames(stratum_areas) <- NULL
  
  # Final areas.
  a1 <- MargAUC(
    areas = data1_areas$area,
    ses = data1_areas$se_area,
    weights = weights
  )
  a0 <- MargAUC(
    areas = data0_areas$area,
    ses = data0_areas$se_area,
    weights = weights
  )
  
  marg_areas <- rbind(a0, a1)
  marg_areas <- cbind(
    'arm' = c(0, 1),
    'n' = c(sum(data0_sizes), sum(data1_sizes)),
    'tau' = tau,
    marg_areas
  )
  
  # Difference and ratio.
  contrasts <- ContrastAreas(
    area1 = marg_areas$area[marg_areas$arm == 1],
    se1 = marg_areas$se[marg_areas$arm == 1],
    area0 = marg_areas$area[marg_areas$arm == 0],
    se0 = marg_areas$se[marg_areas$arm == 0],
    alpha = alpha
  )
  
  # Output
  if(return_areas){
    
    # Average curves.
    aux3 <- function(x){CalcMCF(x$time, x$status, x$idx)}
    data1_curves <- lapply(data1_strata, aux3)
    data0_curves <- lapply(data0_strata, aux3)
    
    data1_mcf <- AvgMCF(data1_curves, weights)
    data0_mcf <- AvgMCF(data0_curves, weights)
    
    data1_mcf$Arm = 1
    data0_mcf$Arm = 0
    
    avg_mcf <- rbind(data0_mcf, data1_mcf)
    
    # Weights.
    weights_df <- data.frame(
      'stratum' = names(data1_sizes),
      'weight' = weights,
      'n' = data0_sizes + data1_sizes,
      'n0' = data0_sizes,
      'n1' = data1_sizes
    )

    # Outputs.
    out <- list(
      'avg_mcf' = avg_mcf,
      'contrasts' = contrasts,
      'marg_areas' = marg_areas,
      'stratum_areas' = stratum_areas,
      'weights' = weights_df
    )
  } else {
    out <- contrasts
  }
  return(out)
}
