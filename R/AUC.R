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


#' Confidence Intervals for the Difference and Ratio of AUCs.
#' 
#' Confidence intervals for the difference and ratio of areas under the mean cumulative
#' count curves, comparing treatment (arm = 1) with reference (arm = 0).
#' 
#' @param time Event time.
#' @param status Status, coded as 1 for the recurrent event, 0 otherwise. 
#' @param arm Arm, coded as 1 for treatment, 0 for reference. 
#' @param idx Subject index. 
#' @param tau Truncation time. 
#' @param B Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom reda mcf Recur
#' @importFrom stats quantile 
#' @export 
#' @return A data.frame containing these columns:
#' \describe{
#'   \item{Time}{Truncation time.}
#'   \item{Arm0}{AUC for arm 0.}
#'   \item{Arm1}{AUC for arm 1.}
#'   \item{Contrast}{'A1-A0' for the difference and 'A1/A0' for the ratio.}
#'   \item{Estimate}{The estimated difference and ratio of AUCs.}
#'   \item{L}{Confidence lower bound.}
#'   \item{U}{Confidence upper bound.}
#' }

AUCCI <- function(time, status, arm, idx, tau, B = 2000, alpha = 0.05){
  
  # Data.
  data <- data.frame(time, status, arm, idx)
  
  # Split data.
  data0 <- data[data$arm == 0, ]
  data1 <- data[data$arm == 1, ]
  
  # Fit mean cumulative function.
  fit <- mcf(Recur(time = time, id = idx, event = status, check = 'none') ~ arm, data = data)
  fit_tab <- fit@MCF
  fit_tab0 <- fit_tab[fit_tab$arm == 0, ]
  fit_tab1 <- fit_tab[fit_tab$arm == 1, ]
  
  # Areas. 
  a0 <- FindAUC(times = fit_tab0$time, values = fit_tab0$MCF, tau = tau)
  a1 <- FindAUC(times = fit_tab1$time, values = fit_tab1$MCF, tau = tau)
  
  # Obseved difference and ratio.
  diff_obs <- a1 - a0
  ratio_obs <- a1 / a0
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0, idx_offset = max(boot1$idx))
    boot <- rbind(boot1, boot0)
    
    # Fit mean cumulative function.
    fitb <- mcf(Recur(time = time, id = idx, event = status, check = 'none') ~ arm, data = boot)
    fitb_tab <- fitb@MCF
    fitb_tab0 <- fitb_tab[fitb_tab$arm == 0, ]
    fitb_tab1 <- fitb_tab[fitb_tab$arm == 1, ]
    
    # Areas. 
    ab0 <- FindAUC(times = fitb_tab0$time, values = fitb_tab0$MCF, tau = tau)
    ab1 <- FindAUC(times = fitb_tab1$time, values = fitb_tab1$MCF, tau = tau)
    
    # Bootstrap difference and ratio.
    diff_boot <- ab1 - ab0
    ratio_boot <- ab1 / ab0
    out <- c(diff_boot, ratio_boot)
    return(out)
  }
  
  boot <- lapply(seq(1:B), aux)
  boot <- do.call(rbind, boot)
  
  # Confidence interval
  alpha2 <- alpha / 2
  L <- apply(boot, 2, function(x){as.numeric(quantile(x, alpha2, na.rm = TRUE))})
  U <- apply(boot, 2, function(x){as.numeric(quantile(x, 1 - alpha2, na.rm = TRUE))})
  
  # Output
  out <- data.frame('Time' = c(tau, tau), 'Arm0' = c(a0, a0), 'Arm1' = c(a1, a1))
  out$Contrast <- c('A1-A0', 'A1/A0')
  out$Estimate <- c(diff_obs, ratio_obs)
  out$L <- L
  out$U <- U
  return(out)
}


#' P-values for the Difference and Ratio of AUCs. 
#' 
#' P-values for the difference and ratio of areas under the mean cumulative
#' count curves, comparing treatment (arm = 1) with control (arm = 0).
#' 
#' @param time Event time.
#' @param status Status, coded as 1 for the recurrent event, 0 otherwise. 
#' @param arm Arm, coded as 1 for treatment, 0 for reference. 
#' @param idx Subject index. 
#' @param tau Truncation time. 
#' @param B Permutations.
#' @param alpha Alpha level.
#' @importFrom reda mcf Recur
#' @export 
#' @return A data.frame containing these columns:
#' \describe{
#'   \item{Time}{Truncation time.}
#'   \item{Arm0}{AUC for arm 0.}
#'   \item{Arm1}{AUC for arm 1.}
#'   \item{Contrast}{'A1-A0' for the difference and 'A1/A0' for the ratio.}
#'   \item{Estimate}{The estimated difference and ratio of AUCs.}
#'   \item{P}{Permutation p-value.}
#' }

AUCP <- function(time, status, arm, idx, tau, B = 2000, alpha = 0.05){
  
  # Data.
  data <- data.frame(time, status, arm, idx)
  
  # Split data.
  data0 <- data[data$arm == 0, ]
  data1 <- data[data$arm == 1, ]
  
  # Fit mean cumulative function.
  fit <- mcf(Recur(time = time, id = idx, event = status, check = 'none') ~ arm, data = data)
  fit_tab <- fit@MCF
  fit_tab0 <- fit_tab[fit_tab$arm == 0, ]
  fit_tab1 <- fit_tab[fit_tab$arm == 1, ]
  
  # Areas. 
  a0 <- FindAUC(times = fit_tab0$time, values = fit_tab0$MCF, tau = tau)
  a1 <- FindAUC(times = fit_tab1$time, values = fit_tab1$MCF, tau = tau)
  
  # Obseved difference and ratio.
  diff_obs <- a1 - a0
  ratio_obs <- a1 / a0
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Permutation data.set.
    data_perm <- PermuteData(data)
    
    # Fit mean cumulative function.
    fitb <- mcf(Recur(time = time, id = idx, event = status, check = 'none') ~ arm, data = data_perm)
    fitb_tab <- fitb@MCF
    fitb_tab0 <- fitb_tab[fitb_tab$arm == 0, ]
    fitb_tab1 <- fitb_tab[fitb_tab$arm == 1, ]
    
    # Areas. 
    ab0 <- FindAUC(times = fitb_tab0$time, values = fitb_tab0$MCF, tau = tau)
    ab1 <- FindAUC(times = fitb_tab1$time, values = fitb_tab1$MCF, tau = tau)
    
    # Permutation contrasts. 
    diff_perm <- ab1 - ab0
    ratio_perm <- ab1 / ab0
    
    # Permutation statistics as or more extreme
    out <- c(NA, NA)
    out[1] <- 1 * (abs(diff_perm) >= abs(diff_obs))
    out[2] <- 1 * (abs(log(ratio_perm)) >= abs(log(ratio_obs)))
    return(out)
  }
  
  perm <- lapply(seq(1:B), aux)
  perm <- do.call(rbind, perm)
  
  # Permutation p-value.
  perm <- rbind(c(1, 1), perm)
  
  # Output
  out <- data.frame('Time' = c(tau, tau), 'Arm0' = c(a0, a0), 'Arm1' = c(a1, a1))
  out$Contrast <- c('A1-A0', 'A1/A0')
  out$Estimate <- c(diff_obs, ratio_obs)
  out$P <- apply(perm, 2, mean)
  return(out)
}