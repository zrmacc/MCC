% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCF.R
\name{CalcVarMCF}
\alias{CalcVarMCF}
\title{Calculate variance of MCF.}
\usage{
CalcVarMCF(mcf, time, status, idx)
}
\arguments{
\item{mcf}{Tabulated MCF as returned by \code{\link{CalcMCF}}.}

\item{time}{Observation times.}

\item{status}{Status indicators.}

\item{idx}{Unique subject index.}
}
\value{
Numeric vector of estimated variances at each distinct observation
  time in `mcf`.
}
\description{
Calculate variance of MCF.
}
