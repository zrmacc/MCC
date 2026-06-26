# Purpose: Data generation.
# Updated: 2022-05-15

#' Resolve death and event frailty variances from legacy alias.
#'
#' @param n_event_types Number of recurrent event types.
#' @param frailty_variance Legacy frailty variance alias.
#' @param death_frailty_variance Explicit death frailty variance or NULL.
#' @param event_frailty_variance Explicit event frailty variance or NULL.
#' @return List with `death_frailty_variance` and `event_frailty_variance`.
#' @noRd
.ResolveFrailtyVariances <- function(
    n_event_types,
    frailty_variance,
    death_frailty_variance,
    event_frailty_variance
) {
  if (n_event_types == 1L) {
    if (!is.null(death_frailty_variance)) {
      death_var <- death_frailty_variance
    } else {
      death_var <- frailty_variance
    }
    if (!is.null(event_frailty_variance) && event_frailty_variance > 0) {
      warning(
        "event_frailty_variance is ignored when n_event_types is 1.",
        call. = FALSE
      )
    }
    return(list(
      death_frailty_variance = death_var,
      event_frailty_variance = 0
    ))
  }

  if (frailty_variance > 0) {
    if (is.null(death_frailty_variance) || is.null(event_frailty_variance)) {
      warning(
        paste(
          "Use death_frailty_variance and event_frailty_variance to specify",
          "death and event frailty variances separately."
        ),
        call. = FALSE
      )
    }
  }

  death_var <- death_frailty_variance
  if (is.null(death_var)) {
    if (frailty_variance > 0) {
      death_var <- frailty_variance
    } else {
      death_var <- 0
    }
  }

  event_var <- event_frailty_variance
  if (is.null(event_var)) {
    if (frailty_variance > 0) {
      event_var <- frailty_variance
    } else {
      event_var <- 0
    }
  }

  return(list(
    death_frailty_variance = death_var,
    event_frailty_variance = event_var
  ))
}

#' Draw mean-one gamma frailties.
#'
#' @param n Number of subjects.
#' @param variance Frailty variance.
#' @return Numeric vector of length `n`.
#' @noRd
.DrawGammaFrailty <- function(n, variance) {
  if (variance > 0) {
    theta <- 1 / variance
    out <- stats::rgamma(n = n, shape = theta, rate = theta)
    return(out)
  }
  return(rep(1, n))
}

#' Build subject-specific baseline recurrent event rates.
#'
#' @param base_event_rate Baseline event rate(s).
#' @param beta_event Covariate log-rate ratios.
#' @param covariates Design matrix.
#' @param n_event_types Number of event types.
#' @param min_event_rate Minimum event rate.
#' @return Numeric matrix with `nrow(covariates)` rows and `n_event_types` cols.
#' @noRd
.ExpandEventRates <- function(
    base_event_rate,
    beta_event,
    covariates,
    n_event_types,
    min_event_rate
) {
  p <- ncol(covariates)
  n <- nrow(covariates)
  base_event_rate <- rep(base_event_rate, length.out = n_event_types)

  if (is.null(beta_event)) {
    beta_event <- rep(0, p)
  }

  rates <- matrix(NA_real_, nrow = n, ncol = n_event_types)

  if (is.matrix(beta_event)) {
    if (nrow(beta_event) != p) {
      stop("beta_event must have one row per covariate column.", call. = FALSE)
    }
    if (ncol(beta_event) != n_event_types) {
      stop(
        "beta_event must have one column per event type when supplied as a matrix.",
        call. = FALSE
      )
    }
    for (k in seq_len(n_event_types)) {
      rates[, k] <- base_event_rate[k] * exp(covariates %*% beta_event[, k])
    }
  } else {
    if (length(beta_event) != p) {
      stop("beta_event must have one coefficient per covariate column.", call. = FALSE)
    }
    linpred <- as.vector(covariates %*% beta_event)
    for (k in seq_len(n_event_types)) {
      rates[, k] <- base_event_rate[k] * exp(linpred)
    }
  }

  rates <- pmax(rates, min_event_rate)
  return(rates)
}

#' Attach type-specific event rates to long-format simulated data.
#'
#' @param data Simulated long-format data.
#' @param event_rates Subject-specific event rates (n x K).
#' @return Data.frame with `event_rate` on recurrent rows.
#' @noRd
.AttachEventRates <- function(data, event_rates) {
  n_event_types <- ncol(event_rates)
  out <- data
  out$event_rate <- NA_real_
  for (k in seq_len(n_event_types)) {
    pick <- out$status == 1L & out$event_type == k
    out$event_rate[pick] <- event_rates[cbind(out$idx[pick], k)]
  }
  return(out)
}

#' Simulate Recurrent Events Data
#'
#' Simulate recurrent events data using exponential censoring, gap, and death
#' times. Status is coded as 0 for censoring, 1 for event, 2 for death.
#' \itemize{
#'   \item The subject-specific death rate is calculated as death_frailty x
#'   base_death_rate x exp(beta_death x covariates).
#'   \item The subject-specific event rate is calculated as death_frailty x
#'   event_frailty x base_event_rate x exp(beta_event x covariates) when
#'   \code{n_event_types > 1}, and death_frailty x base_event_rate x
#'   exp(beta_event x covariates) when \code{n_event_types = 1}.
#' }
#'
#' @param base_death_rate Baseline arrival rate for the terminal event.
#' @param base_event_rate Baseline arrival rate for recurrent events. A scalar
#'   recycled to length \code{n_event_types}, or a length-\code{n_event_types}
#'   numeric vector.
#' @param beta_death Numeric vector of log rate ratios for the death rate.
#' @param beta_event Numeric vector of log rate ratios shared across event types,
#'   or a \code{p x n_event_types} matrix of type-specific log rate ratios.
#' @param censoring_rate Arrival rate for the censoring time. If zero, everyone is
#'   administratively censored at \code{tau} (no random censoring before \code{tau}).
#' @param covariates Numeric design matrix.
#' @param death_frailty_variance Variance of the gamma frailty linking recurrent
#'   and terminal events. When \code{n_event_types = 1}, \code{frailty_variance}
#'   aliases this argument.
#' @param event_frailty_variance Variance of the shared gamma frailty among
#'   recurrent event types. Ignored when \code{n_event_types = 1}.
#' @param frailty_variance Legacy alias for frailty variance. When
#'   \code{n_event_types = 1}, aliases \code{death_frailty_variance}. When
#'   \code{n_event_types > 1} and positive, maps to both frailty variances
#'   unless explicit values are supplied.
#' @param min_death_rate Minimum subject-specific death rate. Must be positive.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n Number of subjects. Overwritten by `nrow(covariates)` if covariates are provided.
#' @param n_event_types Number of recurrent event types.
#' @param tau Truncation time.
#' @return Data.frame. When \code{n_event_types = 1}, columns match the legacy
#'   format (`idx`, `time`, `status`, covariates, `cens_rate`, `death_rate`,
#'   `event_rate`, `frailty`). When \code{n_event_types > 1}, data are in long
#'   format with an `event_type` column on every row, global terminal rows
#'   (`status` 0 or 2) having `event_type = NA`, plus `death_frailty` and
#'   `event_frailty`.
#' @export
GenData <- function(
    base_death_rate = 0.25,
    base_event_rate = 1.0,
    beta_death = NULL,
    beta_event = NULL,
    censoring_rate = 0.25,
    covariates = NULL,
    death_frailty_variance = NULL,
    event_frailty_variance = NULL,
    frailty_variance = 0.0,
    min_death_rate = 0.05,
    min_event_rate = 0.05,
    n = 100,
    n_event_types = 1L,
    tau = 4
) {

  n_event_types <- as.integer(n_event_types)
  if (n_event_types < 1L) {
    stop("n_event_types must be at least 1.", call. = FALSE)
  }

  frailty_vars <- .ResolveFrailtyVariances(
    n_event_types = n_event_types,
    frailty_variance = frailty_variance,
    death_frailty_variance = death_frailty_variance,
    event_frailty_variance = event_frailty_variance
  )

  # Generate subject-specific data frame.
  if (is.null(covariates)) {
    covariates <- data.matrix(rep(1, n))
  } else {
    covariates <- data.matrix(covariates)
    n <- nrow(covariates)
  }
  df <- data.frame(idx = seq_len(n), covariates)

  if (is.null(beta_death)) {
    beta_death <- rep(0, ncol(covariates))
  }
  df$cens_rate <- censoring_rate
  df$death_rate <- base_death_rate * exp(covariates %*% beta_death)
  df$death_rate <- pmax(df$death_rate, min_death_rate)

  df$death_frailty <- .DrawGammaFrailty(n, frailty_vars$death_frailty_variance)
  df$death_rate <- df$death_frailty * df$death_rate

  if (n_event_types == 1L) {
    if (is.null(beta_event)) {
      beta_event <- rep(0, ncol(covariates))
    }
    df$event_rate <- base_event_rate * exp(covariates %*% beta_event)
    df$event_rate <- pmax(df$event_rate, min_event_rate)
    df$frailty <- df$death_frailty
    df$event_rate <- df$frailty * df$event_rate

    data <- SimDataCpp(
      df$cens_rate,
      df$death_rate,
      df$idx,
      df$event_rate,
      tau
    )

    out <- merge(
      x = data,
      y = df[, c("idx", colnames(covariates), "cens_rate", "death_rate", "event_rate", "frailty"), drop = FALSE],
      by = "idx"
    )
    return(out)
  }

  event_rates <- .ExpandEventRates(
    base_event_rate = base_event_rate,
    beta_event = beta_event,
    covariates = covariates,
    n_event_types = n_event_types,
    min_event_rate = min_event_rate
  )

  df$event_frailty <- .DrawGammaFrailty(n, frailty_vars$event_frailty_variance)
  df$frailty <- df$death_frailty * df$event_frailty
  event_rates <- event_rates * df$death_frailty * df$event_frailty

  data <- SimMvDataCpp(
    df$cens_rate,
    df$death_rate,
    df$idx,
    event_rates,
    tau
  )

  data <- .AttachEventRates(data, event_rates)

  out <- merge(
    x = data,
    y = df[, c(
      "idx",
      colnames(covariates),
      "cens_rate",
      "death_rate",
      "death_frailty",
      "event_frailty",
      "frailty"
    ), drop = FALSE],
    by = "idx"
  )
  out <- out[order(out$idx, out$time, out$event_type, na.last = TRUE), ]
  rownames(out) <- NULL
  return(out)
}
