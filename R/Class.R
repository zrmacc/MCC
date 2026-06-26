# -----------------------------------------------------------------------------
# Class for Augmentation Analysis.
# -----------------------------------------------------------------------------

#' Compare Augmented AUCs Object
#'
#' Defines the object class returned by \code{\link{CompareAugAUCs}}.
#'
#' @slot Areas Area under the MCF for each arm.
#' @slot CIs Confidence intervals.
#' @slot MCF Tabulated Mean cumulative function for each treatment arm.
#' @slot Pvals P-values.
#' @slot Reps List of data.frame containing the bootstrap/permutation replicates.
#' @name CompAugAUCs-class
#' @rdname CompAugAUCs-class
#' @exportClass CompAugAUCs

setClass(
  Class = "CompAugAUCs",
  representation = representation(
    Areas = "data.frame",
    CIs = "data.frame",
    MCF = "data.frame",
    Pvals = "data.frame",
    Reps = "list"
  )
)

#' Print Method for Compare Augmented AUCs Object.
#'
#' Print method for objects of class \code{CompAugAUCs}.
#'
#' @param x An object of class \code{CompAugAUCs}.
#' @param ... Unused.
#' @export

print.CompAugAUCs <- function(x, ...) {
  
  disp <- function(y) {
    out <- y
    if (is.numeric(y)) {
      dec_part <- (y %% 1)
      if (max(dec_part, na.rm = TRUE) > 0) {
        out <- signif(y, digits = 3)
      }
    }
    return(out)
  }

  # Areas.
  cat("Marginal Areas:\n")
  areas <- x@Areas
  areas[, ] <- lapply(areas, disp)
  show(areas)
  cat("\n\n")

  # CIs.
  cat("CIs:\n")
  cis <- x@CIs
  cis[, ] <- lapply(cis, disp)
  show(cis)
  cat("\n\n")

  # P-values.
  cat("P-values:\n")
  pvals <- x@Pvals
  pvals[, ] <- lapply(pvals, disp)
  show(pvals)
  cat("\n\n")
}

# -----------------------------------------------------------------------------
# Show Method
# -----------------------------------------------------------------------------

#' Show Method for Compare Augmented AUCs Object
#'
#' @param object An object of class \code{CompAugAUCs}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "CompAugAUCs"),
  definition = function(object) {
    print.CompAugAUCs(x = object)
  }
)


# -----------------------------------------------------------------------------
# Class for Stratified Analysis.
# -----------------------------------------------------------------------------

#' Compare Stratified AUCs Object
#'
#' Defines the object class returned by \code{\link{CompareStratAUCs}}.
#'
#' @slot StratumAreas Stratum-specific areas for each arm.
#' @slot MargAreas Areas for each arm, marginalized across strata.
#' @slot CIs Confidence intervals.
#' @slot MCF Mean cumulative function for each treatment arm, averaged across
#'   strata.
#' @slot Pvals P-values.
#' @slot Reps List of data.frame containing the bootstrap/permutation replicates.
#' @name CompStratAUCs-class
#' @rdname CompStratAUCs-class
#' @exportClass CompStratAUCs

setClass(
  Class = "CompStratAUCs",
  representation = representation(
    StratumAreas = "data.frame",
    MargAreas = "data.frame",
    CIs = "data.frame",
    MCF = "data.frame",
    Pvals = "data.frame",
    Reps = "list"
  )
)

# -----------------------------------------------------------------------------
# Print Method
# -----------------------------------------------------------------------------

#' Print Method for Compare Stratified AUCs Object.
#'
#' Print method for objects of class \code{CompareStratAUCs}.
#'
#' @param x An object of class \code{CompareStratAUCs}.
#' @param ... Unused.
#' @export

print.CompStratAUCs <- function(x, ...) {
  
  disp <- function(y) {
    out <- y
    if (is.numeric(y)) {
      dec_part <- (y %% 1)
      if (max(dec_part, na.rm = TRUE) > 0) {
        out <- signif(y, digits = 3)
      }
      }
    return(out)
  }

  # Areas.
  cat("Marginal Areas:\n")
  areas <- x@MargAreas
  areas[, ] <- lapply(areas, disp)
  show(areas)
  cat("\n\n")

  # CIs.
  cat("CIs:\n")
  cis <- x@CIs
  cis[, ] <- lapply(cis, disp)
  show(cis)
  cat("\n\n")

  # P-values.
  cat("P-values:\n")
  pvals <- x@Pvals
  pvals[, ] <- lapply(pvals, disp)
  show(pvals)
  cat("\n\n")
}

# -----------------------------------------------------------------------------
# Show Method
# -----------------------------------------------------------------------------

#' Show Method for Compare AUCs Object
#'
#' @param object An object of class \code{CompareStratAUCs}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "CompStratAUCs"),
  definition = function(object) {
    print.CompStratAUCs(x = object)
  }
)


# -----------------------------------------------------------------------------
# Class for Multivariate Augmentation Analysis.
# -----------------------------------------------------------------------------

#' Compare Multivariate Augmented AUCs Object
#'
#' Defines the object class returned by \code{\link{CompareMvAUCs}}.
#'
#' @slot Areas Area under the MCF for each arm and event type.
#' @slot CIs Confidence intervals for each event-type contrast.
#' @slot CovDelta Covariance matrix of the contrast vector \eqn{\hat{\Delta}}.
#' @slot CovTheta Named list of arm-specific covariance matrices of
#'   \eqn{\hat{\theta}_a}.
#' @slot MCF Tabulated mean cumulative function for each arm and event type.
#' @slot Pvals P-values for each event-type contrast.
#' @slot Reps List containing bootstrap replicates when requested.
#' @slot Weighted Optional scalar \eqn{w^\top\hat{\Delta}} inference when
#'   contrast \code{process_weights} are supplied to \code{\link{CompareMvAUCs}}.
#' @name CompMvAugAUCs-class
#' @rdname CompMvAugAUCs-class
#' @exportClass CompMvAugAUCs

setClass(
  Class = "CompMvAugAUCs",
  representation = representation(
    Areas = "data.frame",
    CIs = "data.frame",
    CovDelta = "matrix",
    CovTheta = "list",
    MCF = "data.frame",
    Pvals = "data.frame",
    Reps = "list",
    Weighted = "data.frame"
  )
)


#' Print Method for Compare Multivariate Augmented AUCs Object.
#'
#' Print method for objects of class \code{CompMvAugAUCs}.
#'
#' @param x An object of class \code{CompMvAugAUCs}.
#' @param ... Unused.
#' @export

print.CompMvAugAUCs <- function(x, ...) {

  disp <- function(y) {
    out <- y
    if (is.numeric(y)) {
      dec_part <- (y %% 1)
      if (max(dec_part, na.rm = TRUE) > 0) {
        out <- signif(y, digits = 3)
      }
    }
    return(out)
  }

  cat("Marginal Areas:\n")
  areas <- x@Areas
  areas[, ] <- lapply(areas, disp)
  show(areas)
  cat("\n\n")

  cat("CIs:\n")
  cis <- x@CIs
  cis[, ] <- lapply(cis, disp)
  show(cis)
  cat("\n\n")

  cat("P-values:\n")
  pvals <- x@Pvals
  pvals[, ] <- lapply(pvals, disp)
  show(pvals)
  cat("\n\n")

  if (nrow(x@Weighted) > 0) {
    cat("Weighted contrast:\n")
    weighted <- x@Weighted
    weighted[, ] <- lapply(weighted, disp)
    show(weighted)
    cat("\n\n")
  }
}


#' Show Method for Compare Multivariate Augmented AUCs Object
#'
#' @param object An object of class \code{CompMvAugAUCs}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "CompMvAugAUCs"),
  definition = function(object) {
    print.CompMvAugAUCs(x = object)
  }
)
