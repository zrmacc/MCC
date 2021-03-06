% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stratified.R
\name{StratumAUC}
\alias{StratumAUC}
\title{Stratum AUC.}
\usage{
StratumAUC(time, status, idx, tau, calc_var = TRUE)
}
\arguments{
\item{time}{Observation time.}

\item{status}{Status, coded as 0 for censoring, 1 for event, 2 for death.}

\item{idx}{Unique subject index.}

\item{tau}{Truncation time.}

\item{calc_var}{Calculate analytical variance?}
}
\value{
Data.frame containing:
\itemize{
  \item Truncation time 'tau' and estimated 'area'.
  \item Variance 'var_area' and standard error 'se_area', if requested.
}
}
\description{
Calculates the AUC for a single stratum.
}
