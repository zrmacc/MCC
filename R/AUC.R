# Purpose: Functions to perform bootstrap inference on the difference or ratio 
# in AUCs between two mean cumulative count curves. 

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
    y = c(0, values),
    right = TRUE
  )
  area <- integrate(f = g, lower = 0, upper = tau, subdivisions = 2e3)
  return(area$value)
}

# -----------------------------------------------------------------------------

#' Calculate MCF and Find Area
#' 
#' @param time Observation time.
#' @param status Status, coded as 0 for censoring, 1 for event, 2 for death. 
#'   Note that subjects who are neither censored nor die are assumed to
#'   remain at risk throughout the observation period. 
#' @param idx Unique subject index. 
#' @param tau Truncation time. 
#' @return Numeric AUC.

AUC.Area <- function(
  time, 
  status,
  idx,
  tau
) { 
  
  # Fit MCF.
  fit <- CalculateMCF(
    time = time,
    status = status,
    idx = idx
  )
  
  # Find Area.
  area <- FindAUC(
    times = fit$event_times,
    values = fit$mcf,
    tau = tau
  )
  
  # Return area
  return(area)
}

# -----------------------------------------------------------------------------

#' Calculate Test Statistics
#'
#' @param data0 Data.frame containing `time`, `status`, `idx`, `strata` for arm 0.
#' @param data1 Data.frame containing `time`, `status`, `idx`, `strata` for arm 1.
#' @param tau Truncation time.
#' @param return_areas Return the AUCs?
#' @importFrom stats weighted.mean
#' @return If `return_areas = TRUE`, list containing:
#' \itemize{
#'   \item `areas` under the mean cumulative count curve for each arm.
#'   \item `stats`, including the difference and ratio of areas.
#' }
#'  If `return_areas = FALSE`, only `stats` is returned. 

AUC.Stats <- function(data1, data0, tau, return_areas = FALSE) {
  
  # Partition by strata.
  data1_strata <- split(data1, data1$strata, drop = TRUE)
  data0_strata <- split(data0, data0$strata, drop = TRUE)
  
  # Stratum sizes.
  aux1 <- function (x) {length(unique(x$idx))}
  data1_sizes <- sapply(data1_strata, aux1)
  data0_sizes <- sapply(data0_strata, aux1)
  sizes <- data1_sizes + data0_sizes
  weights <- sizes / sum(sizes)
  
  # Stratum areas.
  aux2 <- function (x) {AUC.Area(x$time, x$status, x$idx, tau)}
  data1_areas <- sapply(data1_strata, aux2)
  data0_areas <- sapply(data0_strata, aux2)
  
  # Average curves
  aux3 <- function(x) {CalculateMCF(x$time, x$status, x$idx)}
  data1_curves <- lapply(data1_strata, aux3)
  data0_curves <- lapply(data0_strata, aux3)
  
  data1_mcf <- AverageMCF(data1_curves, weights)
  data0_mcf <- AverageMCF(data0_curves, weights)
  
  data1_mcf$Arm = 1
  data0_mcf$Arm = 0
  
  mcf <- rbind(data0_mcf, data1_mcf)
  
  # Final areas.
  a1 <- mean(data1_areas * weights)
  a0 <- mean(data0_areas * weights)
  areas <- data.frame(
    'N0' = sum(data0_sizes), 
    'Area0' = a0, 
    'N1' = sum(data1_sizes),
    'Area1' = a1
  )
  
  # Difference and ratio
  diff <- a1 - a0
  ratio <- a1 / a0
  stats <- c('diff' = diff, 'ratio' = ratio)
  
  # Output
  if(return_areas){
    weights_df <- data.frame(
      'Stratum' = names(data1_areas),
      'Stratum_weight' = weights,
      'N0' = data0_sizes,
      'Area0' = data0_areas,
      'N1' = data1_sizes,
      'Area1' = data1_areas
    )
    
    out <- list(
      'areas' = areas, 
      'mcf' = mcf,
      'stats' = stats,
      'weights' = weights_df
    )
  } else {
    out <- stats
  }
  return(out)
}