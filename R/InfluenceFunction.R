# Purpose: Influence function at truncation time tau (MCF or AUC).
# Updated: 2026-03-08

#' Calculate Influence Function at Truncation Time
#'
#' Returns the influence function contributions (idx, psi) at the truncation
#' time \code{tau}. One row per subject. The choice of \code{type} determines
#' whether the influence is for the mean cumulative function (MCF) at \code{tau}
#' or for the area under the MCF up to \code{tau} (AUC).
#'
#' @param data Data.frame with columns \code{idx}, \code{status}, \code{time}.
#' @param tau Truncation time.
#' @param type Character. \code{"MCF"} for the MCF influence at \code{tau}, or
#'   \code{"AUC"} for the AUC influence (area under the MCF up to \code{tau}).
#' @param mcf Tabulated MCF as returned by \code{\link{CalcMCF}}. If \code{NULL},
#'   computed from \code{data}.
#' @param weights Optional numeric vector of weights (one per row of \code{data}).
#' @return Data.frame with columns \code{idx}, \code{psi} (one row per subject).
#' @export
InfluenceFunction <- function(
    data,
    tau,
    type = c("MCF", "AUC"),
    mcf = NULL,
    weights = NULL
) {

  type <- match.arg(type)

  if (is.null(weights)) {
    weights <- rep(1, nrow(data))
  }
  data$weights <- weights

  if (is.null(mcf)) {
    mcf <- CalcMCF(
      idx = data$idx,
      status = data$status,
      time = data$time,
      weights = data$weights,
      calc_var = FALSE
    )
  }

  if (type == "AUC") {
    out <- PsiAUC(
      event_rate = mcf$weighted_event_rate,
      grid_time = mcf$time,
      idx = data$idx,
      haz = mcf$haz,
      nar = mcf$nar,
      status = data$status,
      surv = mcf$surv,
      tau = tau,
      time = data$time,
      weights = data$weights
    )
  } else {
    prop_risk <- mcf$nar / mcf$nar[1]
    out <- PsiMCF(
      idx = data$idx,
      event_rate = mcf$weighted_event_rate,
      haz = mcf$haz,
      mcf = mcf$mcf,
      prop_risk = prop_risk,
      status = data$status,
      surv = mcf$surv,
      time = data$time,
      weights = data$weights,
      tau = tau
    )
  }
  return(out)
}
