#' Compare AUCs Object
#'
#' Defines the object class returned by \code{\link{CompareAUCs}}.
#'
#' @slot Areas Overall sample sizes and AUCs for each arm.
#' @slot CIs Equi-tailed and highest density confidence intervals.
#' @slot Curves Mean cumulative function for each treatment arm, averaged across
#'   strata.
#' @slot Pvals Bootstrap and permutation p-values.
#' @slot Reps Bootstrap and permutation statistics. 
#' @slot Weights Per-stratum weights and areas. 
#' @name compAUCs-class
#' @rdname compAUCs-class
#' @exportClass compAUCs

setClass(
  Class = "compAUCs",
  representation = representation(
   Areas = "data.frame",
   CIs = "data.frame",
   Curves = "data.frame",
   Pvals = "data.frame",
   Reps = "matrix",
   Weights = "data.frame"
  )
)

# -----------------------------------------------------------------------------
# Print Method
# -----------------------------------------------------------------------------

#' Print Method for Compre AUCs Object.
#'
#' Print method for objects of class \code{compAUCs}.
#'
#' @param x An object of class \code{compAUCs}.
#' @param ... Unused.
#' @export

print.compAUCs <- function (x, ...) {
  
  # Areas.
  cat('Areas:\n')
  show(x@Areas)
  cat('\n\n')
  
  # CIs.
  cat('CIs:\n')
  show(x@CIs)
  cat('\n\n')
  
  # P-values.
  cat('P-values:\n')
  show(x@Pvals)
  cat('\n\n')

}

# -----------------------------------------------------------------------------
# Show Method
# -----------------------------------------------------------------------------

#' Show Method for Compare AUCs Object
#'
#' @param object An object of class \code{compAUCs}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "compAUCs"),
  definition = function (object) {print.compAUCs(x = object)}
)

