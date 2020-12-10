# -----------------------------------------------------------------------------
# Augmented.
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

print.CompAugAUCs <- function (x, ...) {
  
  disp <- function(y) {
    if (is.numeric(y)) {
      out <- signif(y, digits = 3)
    } else {
      out <- y
    }
    return(out)
  }
  
  # Areas.
  cat('Marginal Areas:\n')
  areas <- x@Areas
  areas[, ] <- lapply(areas, disp)
  show(areas)
  cat('\n\n')
  
  # CIs.
  cat('CIs:\n')
  cis <- x@CIs
  cis[, ] <- lapply(cis, disp)
  show(cis)
  cat('\n\n')
  
  # P-values.
  cat('P-values:\n')
  pvals <- x@Pvals
  pvals[, ] <- lapply(pvals, disp)
  show(pvals)
  cat('\n\n')
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
  definition = function (object) {print.CompAugAUCs(x = object)}
)


# -----------------------------------------------------------------------------
# Stratified.
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

#' Print Method for Compre AUCs Object.
#'
#' Print method for objects of class \code{CompareStratAUCs}.
#'
#' @param x An object of class \code{CompareStratAUCs}.
#' @param ... Unused.
#' @export

print.CompStratAUCs <- function (x, ...) {
  
  disp <- function(y) {
    if (is.numeric(y)) {
      out <- signif(y, digits = 3)
    } else {
      out <- y
    }
    return(out)
  }
  
  # Areas.
  cat('Marginal Areas:\n')
  areas <- x@MargAreas
  areas[, ] <- lapply(areas, disp)
  show(areas)
  cat('\n\n')
  
  # CIs.
  cat('CIs:\n')
  cis <- x@CIs
  cis[, ] <- lapply(cis, disp)
  show(cis)
  cat('\n\n')
  
  # P-values.
  cat('P-values:\n')
  pvals <- x@Pvals
  pvals[, ] <- lapply(pvals, disp)
  show(pvals)
  cat('\n\n')

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
  definition = function (object) {print.CompStratAUCs(x = object)}
)

