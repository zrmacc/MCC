# Purpose: Calculate the average MCF across strata.

#' Calculate Average MCF Curve
#' 
#' Calculates the weighted average of MCF curves for a stratified analysis. 
#' 
#' @param curve_list List of tabulated MCFs as returned by \code{\link{CalcMCF}}.
#' @param weights Numeric vector of weights.
#' @return Data.frame containing `Time` and the averaged MCF `Avg_MCF`.

AvgMCF <- function (curve_list, weights) {
  
  # Extract event times.
  time <- lapply(curve_list, function(x){x$time})
  time <- do.call(c, time)
  time <- sort(unique(time))
  
  # Extract MCFs evaluated on times.
  aux <- function(x){
    g <- stepfun(x$time, c(0, x$mcf), right = FALSE)
    return(g(time))
  }
  mcfs <- lapply(curve_list, aux)
  mcfs <- do.call(cbind, mcfs)
  avg_mcf <- mcfs %*% weights
  
  # Standard errors.
  aux <- function(x){
    g <- stepfun(x$time, c(0, x$var_mcf), right = FALSE)
    return(g(time))
  }
  vars <- lapply(curve_list, aux)
  vars <- do.call(cbind, vars)
  avg_var <- vars %*% (weights^2)
  
  # Output table.
  out <- data.frame(
    'time' = time,
    'mcf' = avg_mcf,
    'var_mcf' = avg_var
  )
  out$se_mcf <- sqrt(out$var_mcf)
  return(out)
}

