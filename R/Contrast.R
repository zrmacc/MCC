# Purpose: Inference on the difference and ratio of AUCs.
# Updated: 2024-02-24

#' Contrast AUCs
#' 
#' Calculate the difference and ratio of areas under the cumulative count curve,
#' comparing two treatment arms.
#' 
#' @param marg_areas Marginal areas, including area and SE for each arm.
#' @param alpha Type I error.
#' @return Data.frame containing the difference and ratio of areas, together
#'   with the confidence interval and p-value.
ContrastAreas <- function(
  marg_areas,
  alpha = 0.05
) {
  
  # Critical value.
  crit <- stats::qnorm(p = 1 - alpha / 2)
  
  # Single-arm case.
  if (nrow(marg_areas) == 1) {
    out <- data.frame(
      contrast = "A0",
      observed = marg_areas$area,
      se = marg_areas$se
    )
    out$lower <- out$observed - crit * out$se
    out$upper <- out$observed + crit * out$se
    out$p <- 2 * stats::pnorm(q = out$observed / out$se, lower.tail = FALSE)
    return(out)
  }
  
  # Unpack.
  area0 <- marg_areas$area[marg_areas$arm == 0]
  area1 <- marg_areas$area[marg_areas$arm == 1]
  
  se0 <- marg_areas$se[marg_areas$arm == 0]
  se1 <- marg_areas$se[marg_areas$arm == 1]
  
  # Difference.
  delta <- area1 - area0
  rho <- area1 / area0
  
  # Output.
  out <- data.frame(
    contrast = c("A1-A0", "A1/A0"),
    observed = c(delta, rho)
  )
  
  if (!is.null(se1) & !is.null(se0)) {
    
    # Inference for delta.
    se_diff <- sqrt(se1^2 + se0^2)
    delta_lower <- delta - crit * se_diff
    delta_upper <- delta + crit * se_diff
    delta_p <- 2 * stats::pnorm(
      q = abs(delta) / se_diff, lower.tail = FALSE)
    
    # Inference for rho.
    rho <- area1 / area0 
    se_rho_log <- sqrt(se1^2 / area1^2 + se0^2 / area0^2)
    rho_lower <- rho * exp(- crit * se_rho_log)
    rho_upper <- rho * exp(+ crit * se_rho_log)
    rho_p <- 2 * stats::pnorm(
      q = abs(log(rho)) / se_rho_log, lower.tail = FALSE)
    
    # Output.
    out$se <- c(se_diff, rho * se_rho_log)
    out$lower <- c(delta_lower, rho_lower)
    out$upper <- c(delta_upper, rho_upper)
    out$p <- c(delta_p, rho_p)
  }
  
  return(out)
}
