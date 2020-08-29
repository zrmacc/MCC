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
  
  # Stratum areas.
  aux2 <- function (x) {AUC.Area(x$time, x$status, x$idx, tau)}
  data1_areas <- sapply(data1_strata, aux2)
  data0_areas <- sapply(data0_strata, aux2)
  
  # Final areas.
  a1 <- weighted.mean(x = data1_areas, w = data1_sizes)
  a0 <- weighted.mean(x = data0_areas, w = data0_sizes)
  areas <- c('a1' = a1, 'a0' = a0)
  
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

# -----------------------------------------------------------------------------

#' Bootstrap Inference on the Area Under the Cumulative Count Curve
#'
#' Confidence intervals and p-values for the difference and ratio of areas under
#' the mean cumulative count curves, comparing treatment (arm = 1) with
#' reference (arm = 0).
#' 
#' Two methods of p-value calculation are available. For '2-sided', treatment 
#' assignments are permuted on each bootstrap iteration, and the p-value is the
#' proportion of the *null* test statistics that are as or more extreme than
#' the *observed* test statistics. For '1-sided', the p-value is the proportion
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
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @param pval_type Either '2-sided' or '1-sided'.
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
#' @examples 
#' \donttest{
#' data(mcc_data)
#' set.seed(100)
#' aucs <- CompareAUCs(
#'   time = mcc_data$time,
#'   status = mcc_data$status,
#'   arm = mcc_data$arm,
#'   idx = mcc_data$idx,
#'   tau = 60,
#'   strata = mcc_data$strata,
#'   reps = 100,
#'   alpha = 0.05,
#'   pval_type = '2-sided'
#' )
#' show(aucs)
#' }

CompareAUCs <- function(
  time, 
  status, 
  arm, 
  idx, 
  tau, 
  strata = NULL,
  reps = 2000, 
  alpha = 0.05,
  pval_type = '2-sided'
) { 
  
  # Create single stratum if no strata are provided. 
  if (is.null(strata)){
    strata <- rep(1, length(time))
  }
  
  # Form data.
  data <- data.frame(time, status, arm, idx, strata)
  data$strata <- factor(data$strata)
  
  # Split data. 
  data0 <- data[data$arm == 0, ]
  data1 <- data[data$arm == 1, ]
  n0 <- nrow(data0)
  
  # Observed test stats.
  obs <- AUC.Stats(
    data0 = data0,
    data1 = data1,
    tau = tau, 
    return_areas = TRUE
  ) 
  obs_areas <- obs$areas
  obs_stats <- obs$stats
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot0 <- StratGroupBoot(data0, idx_offset = 0)
    boot1 <- StratGroupBoot(data1, idx_offset = n0)
    
    # Bootstrap statistics
    boot_stats <- AUC.Stats(
      data0 = boot0,
      data1 = boot1,
      tau = tau
    )
    
    if (pval_type == '1-sided') {
      
      # Indicators of being opposite sign. 
      ind <- as.numeric(sign(boot_stats[1]) != sign(obs_stats[1]))
      ind <- rep(ind, 2)
      
    } else if (pval_type == '2-sided') { 
      
      # Permute data.
      perm <- PermData(rbind(boot0, boot1))
      perm0 <- perm[perm$arm == 0, ]
      perm1 <- perm[perm$arm == 1, ]
      
      # Permutation statistics
      perm_stats <- AUC.Stats(
        data0 = perm0,
        data1 = perm1,
        tau = tau
      )
      
      # Indicators of being more extreme.
      ind1 <- as.numeric(abs(perm_stats[1]) >= abs(obs_stats[1]))
      ind2 <- as.numeric(abs(log(perm_stats[2])) >= abs(log(obs_stats[2])))
      ind <- c(ind1, ind2)
    }
    
    # Results
    out <- c(boot_stats, ind)
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
    FUN = function(x) {as.numeric(quantile(x, 1 - alpha2, na.rm = TRUE))}
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


