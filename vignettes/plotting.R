## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(cache = FALSE, fig.align = "center")
suppressPackageStartupMessages({
  library(dplyr)
  library(MCC)
})

## -----------------------------------------------------------------------------
data <- MCC::GenData(n = 100)

## ----fig.align='center'-------------------------------------------------------
q <- MCC::PlotOneSampleMCF(
  data = data,
  color_lab = "Single Arm"
)
q_nar <- MCC::PlotOneSampleNAR(
  data = data,
  x_breaks = seq(from = 0, to = 4),
  y_lab = "Single Arm"
)
q_main <- cowplot::plot_grid(
  plotlist = list(q, q_nar),
  nrow = 2,
  align = "v",
  axis = "l",
  rel_heights = c(3, 1)
)
show(q_main)

## ----fig.align='center'-------------------------------------------------------
q <- MCC::PlotOneSampleAUMCF(
  data = data,
  color_lab = "Single Arm"
)
q_nar <- MCC::PlotOneSampleNAR(
  data = data,
  x_breaks = seq(from = 0, to = 4),
  y_lab = "Single Arm"
)
q_main <- cowplot::plot_grid(
  plotlist = list(q, q_nar),
  nrow = 2,
  align = "v",
  axis = "l",
  rel_heights = c(3, 1)
)
show(q_main)

## -----------------------------------------------------------------------------
covariates <- data.frame(
  arm = c(rep(1, 100), rep(0, 100)),
  covar = stats::rnorm(200)
)
data <- MCC::GenData(
  beta_event = c(log(0.8), log(1.2)),
  covariates = covariates,
  frailty_variance = 0.2,
  tau = 4
)

## ----fig.align='center'-------------------------------------------------------
q <- MCC::PlotMCFs(
  data = data,
  color_labs = c("Control", "Treatment")
)
q_nar <- MCC::PlotNARs(
  data = data,
  x_breaks = seq(from = 0, to = 4),
  y_labs = c("Control", "Treatment")
)
q_main <- cowplot::plot_grid(
  plotlist = list(q, q_nar),
  nrow = 2,
  align = "v",
  axis = "l",
  rel_heights = c(3, 1)
)
show(q_main)

