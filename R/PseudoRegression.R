# Purpose: Second-order-corrected pseudo-value regression.

#' Pseudo-Value Regression with Second-Order-Corrected Inference
#'
#' Fits an identity-link ordinary least-squares regression to MCF or AUMCF
#' pseudo-values and estimates coefficient covariance using the full
#' first-order Hajek projection, including the second-order contribution.
#'
#' @param formula One-sided regression formula, for example
#'   \code{~ arm + age + arm:age}.
#' @param data Long-format recurrent-event data.
#' @param tau Positive truncation time within the observed time range.
#' @param type \code{"MCF"} for the MCF at \code{tau}, or \code{"AUC"} for
#'   the area under the MCF through \code{tau}.
#' @param covariates Optional subject-level data frame containing the subject
#'   identifier and variables used by \code{formula}. If omitted, these
#'   variables are extracted from \code{data} and must be constant by subject.
#' @param cens_after_last Should subjects without an explicit terminal row be
#'   censored after their last observed event?
#' @param idx_name Name of the subject identifier column.
#' @param status_name Name of the status column: 0 censoring, 1 recurrent event,
#'   and 2 death.
#' @param time_name Name of the observation-time column.
#' @param jump_weights Optional finite recurrent-event jump weights, either one
#'   value or one value per row of \code{data}.
#' @param conf_level Confidence level for normal Wald intervals.
#' @return An object of class \code{PseudoReg}. Its \code{vcov} component and
#'   associated methods use the second-order-corrected covariance. The ordinary
#'   \code{lm} fit is available in \code{fit}.
#' @details The covariance uses the empirical first-order Hajek projection
#'   \code{h11 + h12} with HC0 normalization and normal Wald inference. This
#'   procedure assumes completely independent censoring and the regularity
#'   conditions required for the MCF or AUMCF functional expansion.
#' @examples
#' set.seed(10)
#' design <- cbind(intercept = 1, arm = rep(0:1, each = 10))
#' event_data <- GenData(covariates = design, beta_event = c(0, log(0.8)),
#'                       beta_death = c(0, 0), n = 20, tau = 3)
#' fit <- PseudoReg(~ arm, event_data, tau = 2, type = "MCF")
#' coef(fit)
#' vcov(fit)
#' @export
PseudoReg <- function(
    formula,
    data,
    tau,
    type = c("MCF", "AUC"),
    covariates = NULL,
    cens_after_last = TRUE,
    idx_name = "idx",
    status_name = "status",
    time_name = "time",
    jump_weights = NULL,
    conf_level = 0.95
) {
  call <- match.call()
  type <- match.arg(type)
  .ValidatePseudoFormula(formula)
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  required <- c(idx_name, status_name, time_name)
  if (!all(required %in% names(data))) stop("data is missing an index, status, or time column.", call. = FALSE)
  observed_range <- range(data[[time_name]], na.rm = TRUE)
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau < observed_range[1L] || tau > observed_range[2L]) {
    stop("tau must be a positive finite scalar within the observed time range.", call. = FALSE)
  }
  if (length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.null(jump_weights)) {
    if (!length(jump_weights) %in% c(1L, nrow(data)) || any(!is.finite(jump_weights))) {
      stop("jump_weights must be finite and have length one or nrow(data).", call. = FALSE)
    }
    jump_weights <- rep(jump_weights, length.out = nrow(data))
  }

  prepared <- .PreparePseudo(
    data = data,
    cens_after_last = cens_after_last,
    idx_name = idx_name,
    status_name = status_name,
    tau = tau,
    time_name = time_name,
    type = type,
    jump_weights = jump_weights
  )
  subject <- .PseudoSubjectData(formula, data, covariates, idx_name, prepared$pseudo$idx)
  subject$.pseudo_response <- prepared$pseudo$pseudo
  model_formula <- stats::as.formula(
    call("~", as.name(".pseudo_response"), formula[[2L]]),
    env = environment(formula)
  )
  fit <- stats::lm(model_formula, data = subject, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE)
  design <- fit$x
  if (fit$rank != ncol(design)) stop("The regression model matrix is rank deficient.", call. = FALSE)
  if (nrow(design) <= ncol(design)) stop("The regression requires more subjects than coefficients.", call. = FALSE)

  process_data <- prepared$process_data
  last_time <- tapply(process_data$time, process_data$idx, max)
  if (!any(last_time >= tau)) stop("No subject remains at risk at tau.", call. = FALSE)
  map <- prepared$idx_map[order(prepared$idx_map$idx), , drop = FALSE]
  subject_to_internal <- match(map$orig_idx, subject[[idx_name]])
  design_internal <- design[subject_to_internal, , drop = FALSE]
  h12_internal <- H12MCFCpp(
    idx = process_data$idx,
    status = process_data$status,
    time = process_data$time,
    weights = process_data$jump_weights,
    design = design_internal,
    grid_time = prepared$mcf$time,
    nar = prepared$mcf$nar,
    death = prepared$mcf$death,
    event_weighted = prepared$mcf$weighted_event_rate * prepared$mcf$nar,
    surv = prepared$mcf$surv,
    tau = tau,
    integrate = identical(type, "AUC")
  )
  h12 <- h12_internal[match(subject[[idx_name]], map$orig_idx), , drop = FALSE]
  h11 <- design * as.numeric(stats::residuals(fit))
  score <- h11 + h12
  colnames(h11) <- colnames(h12) <- colnames(score) <- colnames(design)
  rownames(h11) <- rownames(h12) <- rownames(score) <- as.character(subject[[idx_name]])

  n <- nrow(design)
  Ahat <- crossprod(design) / n
  score_centered <- sweep(score, 2L, colMeans(score), FUN = "-")
  Sigmahat <- crossprod(score_centered) / n
  Ainv <- solve(Ahat)
  covariance <- Ainv %*% Sigmahat %*% t(Ainv) / n
  dimnames(covariance) <- list(colnames(design), colnames(design))
  estimate <- stats::coef(fit)
  se <- sqrt(diag(covariance))
  z <- estimate / se
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  coefficient_table <- data.frame(
    term = names(estimate), estimate = unname(estimate), std_error = unname(se),
    z = unname(z), p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    lower = unname(estimate - critical * se), upper = unname(estimate + critical * se),
    row.names = NULL, check.names = FALSE
  )
  pseudo <- prepared$pseudo
  pseudo <- pseudo[match(subject[[idx_name]], pseudo$idx), , drop = FALSE]

  structure(list(
    fit = fit, coefficients = coefficient_table, vcov = covariance,
    pseudo = pseudo, h11 = h11, h12 = h12, score = score,
    Ahat = Ahat, Sigmahat = Sigmahat, tau = tau,
    tau_effective = prepared$tau_effective, type = type, formula = formula,
    conf_level = conf_level, call = call
  ), class = "PseudoReg")
}

.ValidatePseudoFormula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 2L) stop("formula must be one-sided.", call. = FALSE)
  if (any(all.names(formula[[2L]]) == ".")) stop("The formula shorthand '.' is not supported.", call. = FALSE)
  formula_terms <- stats::terms(formula, specials = "offset")
  if (length(attr(formula_terms, "specials")$offset)) stop("Offsets are not supported.", call. = FALSE)
  invisible(TRUE)
}

.PseudoSubjectData <- function(formula, data, covariates, idx_name, subject_ids) {
  variables <- all.vars(formula)
  source <- if (is.null(covariates)) data else covariates
  if (!is.data.frame(source) || !all(c(idx_name, variables) %in% names(source))) {
    stop("The subject identifier and all formula variables must be present in the covariate source.", call. = FALSE)
  }
  if (anyNA(source[, c(idx_name, variables), drop = FALSE])) stop("Missing subject identifiers or model variables are not supported.", call. = FALSE)
  source_ids <- source[[idx_name]]
  if (!setequal(as.character(unique(source_ids)), as.character(subject_ids))) stop("Covariate and event-data subject identifiers do not match.", call. = FALSE)
  if (!is.null(covariates) && anyDuplicated(source_ids)) stop("covariates must contain exactly one row per subject.", call. = FALSE)
  if (is.null(covariates)) {
    for (variable in variables) {
      values <- split(source[[variable]], source_ids)
      if (any(vapply(values, function(x) length(unique(x)) != 1L, logical(1)))) {
        stop(paste0("Covariate '", variable, "' is not constant within subject."), call. = FALSE)
      }
    }
  }
  first <- match(subject_ids, source_ids)
  source[first, c(idx_name, variables), drop = FALSE]
}

#' @export
coef.PseudoReg <- function(object, ...) stats::coef(object$fit)

#' @export
vcov.PseudoReg <- function(object, ...) object$vcov

#' @export
confint.PseudoReg <- function(object, parm, level = object$conf_level, ...) {
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must lie strictly between zero and one.", call. = FALSE)
  }
  estimate <- stats::coef(object$fit)
  se <- sqrt(diag(object$vcov))
  if (missing(parm)) parm <- names(estimate)
  critical <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(estimate - critical * se, estimate + critical * se)
  colnames(out) <- paste0(c((1 - level) / 2, 1 - (1 - level) / 2) * 100, "%")
  out[parm, , drop = FALSE]
}

#' @export
summary.PseudoReg <- function(object, ...) {
  structure(list(call = object$call, type = object$type, tau = object$tau,
                 coefficients = object$coefficients, vcov = object$vcov),
            class = "summary.PseudoReg")
}

#' @export
print.PseudoReg <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.summary.PseudoReg <- function(x, ...) {
  cat("Second-order-corrected pseudo-value regression\n")
  cat("Estimand:", if (x$type == "AUC") "AUMCF" else "MCF", " at tau =", x$tau, "\n\n")
  print(x$coefficients, row.names = FALSE, ...)
  invisible(x)
}
