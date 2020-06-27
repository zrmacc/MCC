#' Sample Repeated Events Data
#'
#' Synthetic time-to event data for two arms of 100 patients, formatted
#' as expected by this package. 
#'
#' @docType data
#' @usage data(mcc_data)
#' @format A data.frame containing four fields:
#' \describe{
#'   \item{idx}{Subject index, 1 through 200.}
#'   \item{time}{Time-to event, integer between 1 and 64.}
#'   \item{status}{Event type, 1 for the recurrent event, 0 otherwise.}
#'   \item{arm}{Treatment arm, 0 for reference, 1 for treatment.}
#' }
"mcc_data"