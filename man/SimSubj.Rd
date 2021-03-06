% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\name{SimSubj}
\alias{SimSubj}
\title{Simulate Data for a Single Subject.}
\usage{
SimSubj(idx, censoring_rate, event_rate, death_rate, tau)
}
\arguments{
\item{idx}{Subject index.}

\item{censoring_rate}{Rate of censoring.}

\item{event_rate}{Rate of events.}

\item{death_rate}{Rate of terminal events.}

\item{tau}{Truncation time.}
}
\value{
Recurrent event data for a single subject.
}
\description{
Simulate Data for a Single Subject.
}
