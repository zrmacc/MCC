% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Resampling.R
\name{StratGroupBoot}
\alias{StratGroupBoot}
\title{Stratified, Grouped Bootstrap.}
\usage{
StratGroupBoot(data)
}
\arguments{
\item{data}{Data.frame containing `idx` and `strata`.}
}
\value{
Bootstrapped data.frame.
}
\description{
Perform grouped bootstrap within levels of a stratification factor.
}
