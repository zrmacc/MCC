% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{SingleArmAUC}
\alias{SingleArmAUC}
\title{Single Arm Area Under the Cumulative Count Curve}
\usage{
SingleArmAUC(
  data,
  alpha = 0.05,
  boot = FALSE,
  cens_after_last = TRUE,
  idx_name = "idx",
  reps = 2000,
  status_name = "status",
  strata = NULL,
  tau = NULL,
  time_name = "time",
  weights = NULL
)
}
\arguments{
\item{data}{Data.frame.}

\item{alpha}{Type I error level.}

\item{boot}{Logical, construct bootstrap confidence intervals (slow)?}

\item{cens_after_last}{Should subjects who lack an explicit censoring time
be censored after their last observed event?}

\item{idx_name}{Name of column containing a unique subject index.}

\item{reps}{Number of replicates for bootstrap inference.}

\item{status_name}{Name of column containing the status. Must be coded as 0 for
censoring, 1 for event, 2 for death. Each subject should have an 
observation-terminating event, either censoring or death.}

\item{strata}{Optional stratification factor.}

\item{tau}{Numeric truncation time.}

\item{time_name}{Name of column containing the observation time.}

\item{weights}{Optional column of weights, controlling the size of the jump
in the cumulative count curve at times with status == 1.}
}
\value{
Object of class compAUCs with these slots:
\itemize{
  \item `@Areas`: The AUC for each arm.
  \item `@CIs`: Observed difference and ratio in areas with confidence intervals.
  \item `@Curves`: Mean cumulative count curve for each arm; averaged across strata
    if present.
  \item `@Pvals`: Bootstrap and permutation p-values.
  \item `@Reps`: Bootstrap and permutation realizations of the test statistics.
  \item `@Weights`: Per-stratum weights and AUCs.
}
}
\description{
Single Arm Area Under the Cumulative Count Curve
}
\examples{
# Simulate data set.
covar <- data.frame(strata = rep(c(1, 2), each = 50))
data <- GenData(beta_event = log(0.5), covariates = covar)
# Calculate AUC.
auc <- SingleArmAUC(data, strata = data$strata, tau = 2)
\donttest{auc <- SingleArmAUC(data, boot = TRUE, reps = 100, strata = data$strata, tau = 2)}
show(auc)
}
