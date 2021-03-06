% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Contrast.R
\name{ContrastAreas}
\alias{ContrastAreas}
\title{Contrast AUCs}
\usage{
ContrastAreas(marg_areas, alpha)
}
\arguments{
\item{marg_areas}{Marginal areas, including area and SE for each arm.}

\item{alpha}{Type I error.}
}
\value{
Data.frame containing:
\itemize{
  \item 'Contrast' and estimate 'Est'.
  \item Lower 'L' and upper 'U' confidence bounds.
  \item 'P' value.
}
}
\description{
Calculate the difference and ratio of areas under the cumulative count curve,
comparing two treatment arms.
}
