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


#' Calculate Test Statistics
#'
#' @param data Data.frame containing `time`, `status`, `arm`, 
#'   and subject `idx`. 
#' @param tau Truncation time.
#' @param return_areas Return the AUCs?
#' @importFrom reda mcf Recur
#' @return If `return_areas = TRUE`, list containing:
#' \itemize{
#'   \item `areas` under the mean cumulative count curve for each arm.
#'   \item `stats`, including the difference and ratio of areas.
#' }
#'  If `return_areas = FALSE`, only `stats` is returned. 

AUC.Stats <- function(data, tau, return_areas = FALSE) {
  
  # Fit mean cumulative function.
  fit <- mcf(Recur(time = time, id = idx, event = status, check = 'none') ~ arm, data = data)
  fit_tab <- fit@MCF
  fit_tab0 <- fit_tab[fit_tab$arm == 0, ]
  fit_tab1 <- fit_tab[fit_tab$arm == 1, ]
  
  # Areas. 
  a0 <- FindAUC(times = fit_tab0$time, values = fit_tab0$MCF, tau = tau)
  a1 <- FindAUC(times = fit_tab1$time, values = fit_tab1$MCF, tau = tau)
  areas <- c('a0' = a0, 'a1' = a1)
  
  # Difference and ratio
  diff <- a1 - a0
  ratio <- a1 / a0
  stats <- c('diff' = diff, 'ratio' = ratio)
  
  # Output
  if(return_areas){
    out <- list('areas' = areas, 'stats' = stats)
  } else {
    out <- stats
  }
  return(out)
}


#' Bootstrap Inference on the Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0).
#' 
#' @param time Event time.
#' @param status Status, coded as 1 for the recurrent event, 0 otherwise. 
#' @param arm Arm, coded as 1 for treatment, 0 for reference. 
#' @param idx Subject index. 
#' @param tau Truncation time. 
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom stats quantile 
#' @export 
#' @return Data.frame containing these columns:
#' \describe{
#'   \item{Time}{Truncation time.}
#'   \item{Arm0}{AUC for arm 0.}
#'   \item{Arm1}{AUC for arm 1.}
#'   \item{Contrast}{'A1-A0' for the difference and 'A1/A0' for the ratio.}
#'   \item{Estimate}{The estimated difference and ratio of AUCs.}
#'   \item{L}{Confidence lower bound.}
#'   \item{U}{Confidence upper bound.}
#'   \item{P}{P-value.}
#' }

CompareAUCs <- function(
  time, 
  status, 
  arm, 
  idx, 
  tau, 
  reps = 2000, 
  alpha = 0.05
) { 
  
  # Data.
  data <- data.frame(time, status, arm, idx)
  
  # Observed test stats.
  obs <- AUC.Stats(data = data, tau = tau, return_areas = TRUE) 
  obs_areas <- obs$areas
  obs_stats <- obs$stats
  
  # Split data.
  data0 <- data[data$arm == 0, ]
  data1 <- data[data$arm == 1, ]
  n1 <- nrow(data1)
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0, idx_offset = n1)
    boot <- rbind(boot1, boot0)
    
    # Bootstrap statistics
    boot_stats <- AUC.Stats(data = boot, tau = tau)
    
    # Permute data.
    perm <- PermData(boot)
    
    # Permutation statistics
    perm_stats <- AUC.Stats(data = perm, tau = tau)
    
    # Results
    out <- boot_stats
    out[3] <- 1 * abs(perm_stats[1]) >= abs(obs_stats[1])
    out[4] <- 1 * abs(log(perm_stats[2])) >= abs(log(obs_stats[2]))
    return(out)
  }
  
  sim <- lapply(seq(1:reps), aux)
  sim <- do.call(rbind, sim)
  
  # Confidence interval
  alpha2 <- alpha / 2
  lower <- apply(
    X = sim[, 1:2], 
    MARGIN = 2, 
    FUN = function(x) {as.numeric(quantile(x, alpha2, na.rm = TRUE))}
  )
  upper <- apply(
    X = sim[, 1:2], 
    MARGIN = 2, 
    FUN = function(x){as.numeric(quantile(x, 1 - alpha2, na.rm = TRUE))}
  )
  
  # P-values
  pval <- apply(
    X = rbind(sim[, 3:4], c(1, 1)),
    MARGIN = 2,
    FUN = function(x){as.numeric(mean(x, na.rm = TRUE))}
  )
  
  # Output
  out <- data.frame(
    'Time' = c(tau, tau), 
    'Arm0' = rep(x = obs_areas[1], times = 2),
    'Arm1' = rep(x = obs_areas[2], times = 2)
  )
  out$Contrast <- c('A1-A0', 'A1/A0')
  out$Estimate <- obs_stats
  out$L <- lower
  out$U <- upper
  out$P <- pval
  return(out)
}


