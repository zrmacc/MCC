# Tests for multivariate MCF plotting functions.

suppressPackageStartupMessages({
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  return(x)
}

#' Build simple multivariate long-format data for plotting tests.
#'
#' @param specs List of per-subject specifications.
#' @return Long-format data.frame.
.MvPlotTestData <- function(specs) {
  rows <- lapply(specs, function(s) {
    arm <- s$arm
    idx <- s$idx
    out <- NULL
    for (ev in s$events) {
      out <- rbind(
        out,
        data.frame(
          arm = arm,
          idx = idx,
          time = ev$time,
          status = ev$status,
          event_type = ev$event_type,
          stringsAsFactors = FALSE
        )
      )
    }
    return(out)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}


#' Expect a ggplot or cowplot gtable plot object.
#'
#' @param p Plot object.
#' @noRd
.ExpectMvPlot <- function(p) {
  is_plot <- inherits(p, "ggplot") || inherits(p, "gtable")
  testthat::expect_true(is_plot)
}


test_that("PlotMvMCFs runs on multivariate GenData output.", {
  withr::local_seed(201)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  p <- PlotMvMCFs(data = data, tau = 4)
  .ExpectMvPlot(p)
  expect_error(show(p), NA)
})

test_that("PlotMvMCFs default includes NAR per event type; show_nar=FALSE is MCF-only.", {
  withr::local_seed(208)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  p_combined <- PlotMvMCFs(data = data, tau = 4, x_breaks = 0:4)
  .ExpectMvPlot(p_combined)
  expect_s3_class(p_combined, "ggplot")
  p_mcf <- PlotMvMCFs(data = data, tau = 4, show_nar = FALSE)
  expect_s3_class(p_mcf, "ggplot")
  expect_equal(p_mcf$theme$legend.position, "top")
  expect_s3_class(ggplot2::ggplot_build(p_mcf)$layout$facet, "FacetWrap")
  expect_false(inherits(ggplot2::ggplot_build(p_combined)$layout$facet, "FacetWrap"))
})

test_that("PlotMvMCFs one-arm via which_arm includes NAR per event type.", {
  withr::local_seed(209)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  p_combined <- PlotMvMCFs(
    data = data,
    which_arm = 0,
    tau = 4,
    x_breaks = 0:4,
    color_lab = "Control"
  )
  .ExpectMvPlot(p_combined)
  expect_s3_class(p_combined, "ggplot")
  p_mcf <- PlotMvMCFs(
    data = data,
    which_arm = 0,
    tau = 4,
    show_nar = FALSE,
    color_lab = "Control"
  )
  expect_s3_class(p_mcf, "ggplot")
  expect_equal(p_mcf$theme$legend.position, "top")
  expect_s3_class(ggplot2::ggplot_build(p_mcf)$layout$facet, "FacetWrap")
  expect_false(inherits(ggplot2::ggplot_build(p_combined)$layout$facet, "FacetWrap"))
  expect_error(show(p_combined), NA)
})

test_that("PlotMvMCFs accepts custom event_type_labs.", {
  withr::local_seed(210)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  labs <- c("Hospitalization", "Heart failure")
  p <- PlotMvMCFs(
    data = data,
    tau = 4,
    x_breaks = 0:4,
    event_type_labs = labs,
    show_nar = FALSE
  )
  built <- ggplot2::ggplot_build(p)
  expect_equal(
    as.character(built$layout$layout$event_type_label),
    labs
  )
})

test_that("PlotMvMCFs delegates to PlotMCFs when event_type is absent.", {
  withr::local_seed(202)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(covariates = covariates, tau = 4)
  p_mv <- PlotMvMCFs(data = data, tau = 4, show_nar = FALSE)
  p_uni <- PlotMCFs(data = data, tau = 4, show_nar = FALSE)
  expect_s3_class(p_mv, "ggplot")
  expect_s3_class(p_uni, "ggplot")
})

test_that("PlotMvMCFs show_auc runs for one- and two-arm multivariate data.", {
  withr::local_seed(203)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  expect_error(show(PlotMvMCFs(
    data = data,
    tau = 4,
    show_nar = FALSE,
    show_auc = TRUE
  )), NA)
  cov_one <- data.frame(arm = rep(0, 50))
  data_one <- GenData(
    covariates = cov_one,
    n_event_types = 2L,
    tau = 4
  )
  expect_error(show(PlotMvMCFs(
    data = data_one,
    tau = 4,
    show_nar = FALSE,
    show_auc = TRUE,
    color_lab = "Control"
  )), NA)
})

test_that("PlotMvMCFs show_auc with show_nar runs without error.", {
  withr::local_seed(211)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  p <- PlotMvMCFs(
    data = data,
    tau = 4,
    x_breaks = 0:4,
    show_nar = TRUE,
    show_auc = TRUE
  )
  .ExpectMvPlot(p)
  expect_error(show(p), NA)
})

test_that("PlotMvMCFs supports K = 3 arms with custom colors and labels.", {
  withr::local_seed(212)
  data <- .MvPlotTestData(list(
    list(
      arm = 0, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 1, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 1, event_type = 2L),
        list(time = 4, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 2, idx = 3,
      events = list(
        list(time = 2, status = 1, event_type = 2L),
        list(time = 4, status = 0, event_type = NA)
      )
    )
  ))
  p <- PlotMvMCFs(
    data = data,
    tau = 4,
    x_breaks = 0:4,
    color_labs = c("Placebo", "Low dose", "High dose"),
    colors = c("#C65842", "#6385B8", "#5BA85B")
  )
  .ExpectMvPlot(p)
  expect_error(show(p), NA)
})

test_that("PlotMvMCFs errors when which_arm is not found.", {
  withr::local_seed(213)
  covariates <- data.frame(arm = rep(c(0, 1), each = 25))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  expect_error(
    PlotMvMCFs(data = data, which_arm = 9, tau = 4),
    "not found"
  )
})

test_that("MvMCFCurve and MvNARCurve return one curve per event type.", {
  withr::local_seed(205)
  covariates <- data.frame(arm = rep(0, 30))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  mcf_curves <- MvMCFCurve(data = data)
  nar_curves <- MvNARCurve(data = data)
  expect_length(mcf_curves, 2)
  expect_length(nar_curves, 2)
  expect_true(all(vapply(mcf_curves, inherits, logical(1), "stepfun")))
  expect_true(all(vapply(nar_curves, inherits, logical(1), "stepfun")))
})

test_that("MvMCFPlotFrame returns rows for each event type.", {
  withr::local_seed(206)
  covariates <- data.frame(arm = rep(c(0, 1), each = 20))
  data <- GenData(
    covariates = covariates,
    n_event_types = 2L,
    tau = 4
  )
  norm_mv <- getFromNamespace(".NormMvDataForPlot", "MCC")
  mv_prep <- norm_mv(
    data = data,
    arm_name = "arm",
    cens_after_last = TRUE,
    event_type_name = "event_type",
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
  )
  frame <- MvMCFPlotFrame(mv_prep = mv_prep, tau = 4)
  expect_equal(sort(unique(frame$event_type)), c(1L, 2L))
  expect_true("event_type_label" %in% names(frame))
})

test_that("Typed terminal events produce plots without error.", {
  tau <- 5
  data <- .MvPlotTestData(list(
    list(
      arm = 1, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = 1L),
        list(time = 3, status = 1, event_type = 2L),
        list(time = 5, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 0, idx = 2,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 2, status = 0, event_type = NA)
      )
    )
  ))
  expect_error(PlotMvMCFs(data = data, tau = tau, cens_after_last = FALSE), NA)
})

test_that("At-risk subsetting produces plots without error.", {
  tau <- 5
  data <- .MvPlotTestData(list(
    list(
      arm = 0, idx = 1,
      events = list(
        list(time = 1, status = 1, event_type = 1L),
        list(time = 3, status = 0, event_type = NA)
      )
    ),
    list(
      arm = 1, idx = 2,
      events = list(
        list(time = 2, status = 1, event_type = 1L),
        list(time = 2, status = 1, event_type = 2L),
        list(time = 4, status = 0, event_type = NA)
      )
    )
  ))
  p <- PlotMvMCFs(
    data = data,
    tau = tau,
    cens_after_last = FALSE,
    x_breaks = 0:tau
  )
  .ExpectMvPlot(p)
})
