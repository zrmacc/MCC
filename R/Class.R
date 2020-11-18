#' Compare AUCs Object
#'
#' Defines the object class returned by \code{\link{CompareAUCs}}.
#'
#' @slot StratumAreas Stratum-specific areas for each arm.
#' @slot MargAreas Areas for each arm, marginalized across strata.
#' @slot CIs Confidence intervals.
#' @slot MCF Mean cumulative function for each treatment arm, averaged across
#'   strata.
#' @slot Pvals P-values.
#' @slot Reps List of data.frame containing the bootstrap/permutation replicates.
#' @slot Weights Per-stratum weights and areas. 
#' @name compAUCs-class
#' @rdname compAUCs-class
#' @exportClass compAUCs

setClass(
  Class = "compAUCs",
  representation = representation(
   StratumAreas = "data.frame",
   MargAreas = "data.frame",
   CIs = "data.frame",
   MCF = "data.frame",
   Pvals = "data.frame",
   Reps = "list",
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
#' @param object An object of class \code{compAUCs}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "compAUCs"),
  definition = function (object) {print.compAUCs(x = object)}
)

