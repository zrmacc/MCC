% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCF.R
\name{CalcMCF}
\alias{CalcMCF}
\title{Calculate Mean Cumulative Function}
\usage{
CalcMCF(time, status, idx, calc_var = TRUE)
}
\arguments{
\item{time}{Observation time.}

\item{status}{Status, coded as 0 for censoring, 1 for event, 2 for death.}

\item{idx}{Unique subject index.}

\item{calc_var}{Calculate variance of the MCF?}
}
\value{
Data.frame with these columns:
\itemize{
  \item `times`, distinct observation times.
  \item `censor`, number of censorings.
  \item `death`, number of deaths.
  \item `event`, number of events.
  \item `haz`, instantaneous hazard (of death).
  \item `surv`, survival probability.
  \item `event_rate`, instantaneous event rate.
  \item `mcf`, mean cumulative function. 
}
}
\description{
Tabulates the mean cumulative function. See equation 2.1 of 
<doi:10.1111/j.0006-341X.2000.00554.x>.
}
