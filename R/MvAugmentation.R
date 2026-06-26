# Purpose: Multivariate augmentation components.
# Updated: 2026-06-16

#' Build Stacked Centered Covariate Matrix
#'
#' @param covars Subject-level covariate matrix (n x p).
#' @param xbar_by_k List of length K with p-vectors of covariate means.
#' @return Matrix (n x Kp).
#' @noRd
.MvStackResid <- function(
    covars,
    xbar_by_k
) {

  K <- length(xbar_by_k)
  p <- ncol(covars)
  n <- nrow(covars)
  out <- matrix(0, nrow = n, ncol = K * p)
  for (k in seq_len(K)) {
    cols <- ((k - 1L) * p + 1L):(k * p)
    out[, cols] <- sweep(covars, 2, xbar_by_k[[k]], "-")
  }
  return(out)
}


#' Calculate Multivariate Augmentation Components for One Arm
#'
#' @param covars Subject-level covariate matrix (n x p).
#' @param psi_mat Influence matrix (n x K).
#' @param xbar_by_k Named list of covariate means per event type.
#' @param n_a Arm sample size.
#' @return List with sigma (Kp x Kp) and gamma (Kp x K).
#' @noRd
CalcMvAugCompArm <- function(
    covars,
    psi_mat,
    xbar_by_k,
    n_a
) {

  K <- ncol(psi_mat)
  p <- ncol(covars)
  resid <- .MvStackResid(covars = covars, xbar_by_k = xbar_by_k)
  sigma <- crossprod(resid) / (n_a * n_a)
  gamma <- matrix(0, nrow = K * p, ncol = K)
  for (k in seq_len(K)) {
    gamma[, k] <- crossprod(resid, psi_mat[, k]) / (n_a * n_a)
  }
  out <- list(
    sigma = sigma,
    gamma = gamma,
    xbar_by_k = xbar_by_k
  )
  return(out)
}
