# Purpose: Calculate the average MCF across strata.

#' Calculate Average MCF Curve
#'
#' Calculates the weighted average of MCF curves for a stratified analysis.
#'
#' @param curve_list List of tabulated MCFs as returned by \code{\link{CalcMCF}}.
#' @param weights Numeric vector of weights.
#' @return Data.frame containing `Time` and the averaged MCF `Avg_MCF`.

AvgMCF <- function(curve_list, weights) {

  # Extract event times.
  time <- lapply(curve_list, function(x) {x$time})
  time <- do.call(c, time)
  time <- sort(unique(time))

  # Extract MCFs evaluated on times.
  aux <- function(x) {
    g <- stepfun(x$time, c(0, x$mcf), right = FALSE)
    return(g(time))
  }
  mcfs <- lapply(curve_list, aux)
  mcfs <- do.call(cbind, mcfs)
  avg_mcf <- mcfs %*% weights

  # Standard errors.
  aux <- function(x) {
    g <- stepfun(x$time, c(0, x$var_mcf), right = FALSE)
    return(g(time))
  }
  vars <- lapply(curve_list, aux)
  vars <- do.call(cbind, vars)
  avg_var <- vars %*% (weights^2)

  # Output table.
  out <- data.frame(
    "time" = time,
    "mcf" = avg_mcf,
    "var_mcf" = avg_var
  )
  out$se_mcf <- sqrt(out$var_mcf)
  return(out)
}


# -----------------------------------------------------------------------------
# Calculate Marginal MCF
# -----------------------------------------------------------------------------

#' Calculate Marginal MCF
#' 
#' Calculates the marginal MCF, averaged across strata, with stratum
#' weights proportional to the total number of subjects (across arms)
#' belonging to that stratum. 
#'
#' @param data Data.frame containing: idx, time, status, arm, strata.
#' @importFrom dplyr "%>%"
#' @return Data.frame.
#' @export 

CalcMargMCF <- function(data) {
  
  # Stratum sizes.
  idx <- time <- status <- arm <- strata <- NULL
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
  
  # Marginal MCF for arm 1.
  mcf1 <- data %>%
    dplyr::filter(arm == 1) %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      CalcMCF(time, status, idx),
      .groups = "keep"
    ) %>%
    dplyr::group_split()
  avg_mcf1 <- AvgMCF(mcf1, weights = stratum_sizes$w1[stratum_sizes$w1 != 0])
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
  avg_mcf0 <- AvgMCF(mcf0, weights = stratum_sizes$w0[stratum_sizes$w0 != 0])
  avg_mcf0$arm <- 0
  
  avg_mcf <- rbind(avg_mcf1, avg_mcf0)
  return(avg_mcf)
}
