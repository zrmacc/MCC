#!/usr/bin/env Rscript

# Run after installing the current MCC source tree.
suppressPackageStartupMessages({
  library(MCC)
  library(testthat)
})

test_file("tests/testthat/test-pseudo-regression.R", stop_on_failure = TRUE)

set.seed(20260903)
max_centering_error <- 0
for (type in c("MCF", "AUC")) {
  for (replicate in seq_len(25)) {
    design <- cbind(intercept = 1, A = rep(0:1, each = 20), X = rnorm(40))
    data <- GenData(
      covariates = design,
      beta_event = c(0, .2, -.1),
      beta_death = c(0, 0, 0),
      n = 40,
      tau = 3
    )
    fit <- PseudoReg(~ A + X + A:X, data, tau = 2, type = type)
    stopifnot(all(is.finite(fit$vcov)))
    max_centering_error <- max(max_centering_error, abs(colMeans(fit$h12)))
  }
}
stopifnot(max_centering_error < 1e-10)
message("Validation complete; maximum h12 centering error: ", format(max_centering_error, scientific = TRUE))
