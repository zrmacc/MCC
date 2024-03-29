% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AUC.R
\name{VarAUC}
\alias{VarAUC}
\title{Calculate Variance of AUC}
\usage{
VarAUC(data, tau, mcf = NULL, return_psi = FALSE, weights = NULL)
}
\arguments{
\item{data}{Data.frame containing {idx, status, time}.}

\item{tau}{Truncation time.}

\item{mcf}{Tabulated MCF, if already computed.}

\item{return_psi}{Return influence function contributions?}

\item{weights}{Optional column of weights, controlling the size of the jump
in the cumulative count curve at times with status == 1.}
}
\value{
Numeric variance.
}
\description{
Calculate Variance of AUC
}
