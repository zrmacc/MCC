# Purpose: Univariate mean cumulative function plotting.
# Updated: 2026-06-17

# ------------------------------------------------------------------------------


#' Plot Mean Cumulative Functions
#'
#' Plot mean cumulative functions for one or two treatment arms, with optional
#' numbers at risk below the MCF panel (\code{show_nar = TRUE}) and optional
#' area shading (\code{show_auc = TRUE}).
#'
#' @section Migration:
#' Prior releases returned an MCF-only \code{ggplot} from \code{PlotMCFs()}.
#' The default is now a combined MCF + NAR layout. Pass \code{show_nar = FALSE}
#' for MCF-only output. Legacy functions \code{PlotOneSampleMCF},
#' \code{PlotOneSampleAUMCF}, \code{PlotOneSampleNAR}, \code{PlotAUMCF}, and
#' \code{PlotNARs} remain as deprecated wrappers.
#'
#' @param data Data.frame.
#' @param which_arm Optional arm code to plot a single arm; \code{NULL} plots
#'   all arms present.
#' @param arm_name Name of arm column in data.
#' @param color_labs Legend labels for two-arm plots.
#' @param colors Line and fill colors for two-arm plots.
#' @param ctrl_color Deprecated control color (used when \code{colors} is
#'   \code{NULL}).
#' @param trt_color Deprecated treatment color (used when \code{colors} is
#'   \code{NULL}).
#' @param color Line color when a single arm is plotted.
#' @param color_lab Legend label when a single arm is plotted.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param strata_name Name of stratum column in data.
#' @param tau Truncation time.
#' @param time_name Name of time column in data.
#' @param title Plot title.
#' @param jump_weights Optional column of jump weights for one-arm MCF jumps.
#' @param show_nar Include numbers at risk below the MCF panel?
#' @param show_auc Shade the area under each MCF curve from 0 to \code{tau}?
#' @param show_mcf Include the MCF panel? Set \code{FALSE} for NAR-only output.
#' @param x_breaks X-axis breaks for the NAR row.
#' @param x_labs X-axis tick labels for the NAR row.
#' @param x_max X-axis upper limit for the NAR row.
#' @param x_lim X-axis limits.
#' @param x_name X-axis label.
#' @param y_breaks Y-axis breaks.
#' @param y_lim Y-axis limits.
#' @param y_name Y-axis label.
#' @return A \code{ggplot} object, or a cowplot layout when \code{show_nar =
#'   TRUE} and \code{show_mcf = TRUE}.
#' @importFrom dplyr "%>%"
#' @export
PlotMCFs <- function(
    data,
    which_arm = NULL,
    arm_name = "arm",
    color_labs = c("Ctrl", "Trt"),
    colors = NULL,
    ctrl_color = "#C65842",
    trt_color = "#6385B8",
    color = "#C65842",
    color_lab = "Arm",
    idx_name = "idx",
    status_name = "status",
    strata_name = NULL,
    tau = NULL,
    time_name = "time",
    title = NULL,
    jump_weights = NULL,
    show_nar = TRUE,
    show_auc = FALSE,
    show_mcf = TRUE,
    x_breaks = NULL,
    x_labs = NULL,
    x_max = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  if (!show_mcf && !show_nar) {
    stop("At least one of show_mcf or show_nar must be TRUE.", call. = FALSE)
  }

  data <- .MvSubsetWhichArm(
    data = data,
    which_arm = which_arm,
    arm_name = arm_name
  )
  one_arm <- .UniIsOneArm(data = data, arm_name = arm_name)
  arm_style <- NULL

  if (one_arm) {
    data_norm <- .NormDataForPlot(
      data = data,
      arm_name = NULL,
      strata_name = NULL,
      idx_name = idx_name,
      status_name = status_name,
      time_name = time_name,
      jump_weights = jump_weights
    )
  } else {
    data_norm <- .NormDataForPlot(
      data = data,
      arm_name = arm_name,
      strata_name = strata_name,
      idx_name = idx_name,
      status_name = status_name,
      time_name = time_name,
      jump_weights = jump_weights
    )
    arms <- sort(unique(data_norm[["arm"]]))
    if (!identical(as.integer(arms), c(0L, 1L))) {
      stop(
        "Two-arm plots require exactly arms 0 and 1.",
        call. = FALSE
      )
    }
    arm_style <- .UniResolveArmStyle(
      color_labs = color_labs,
      colors = colors,
      ctrl_color = ctrl_color,
      trt_color = trt_color
    )
  }

  if (is.null(x_lim) || is.null(x_lim[2])) {
    x_max_data <- max(data_norm[["time"]])
  } else {
    x_max_data <- x_lim[2]
  }
  if (is.null(x_max)) {
    x_max <- x_max_data
  }
  if (is.null(tau)) {
    tau <- x_max_data
  }
  if (show_nar && is.null(x_breaks)) {
    x_breaks <- .MvDefaultXBreaks(tau)
  }

  if (show_mcf && show_nar) {
    out <- .UniBuildCombinedPlot(
      one_arm = one_arm,
      data_norm = data_norm,
      arm_style = arm_style,
      tau = tau,
      x_breaks = x_breaks,
      x_labs = x_labs,
      x_max = x_max,
      color = color,
      color_lab = color_lab,
      show_auc = show_auc,
      title = title,
      x_lim = x_lim,
      x_name = x_name,
      y_breaks = y_breaks,
      y_lim = y_lim,
      y_name = y_name
    )
    return(out)
  }

  if (show_mcf && !show_nar) {
    if (one_arm) {
      mcf_df <- .UniOneArmMCFFrame(data = data_norm, tau = tau)
      out <- .UniOneArmMCFPanel(
        mcf_df = mcf_df,
        color = color,
        color_lab = color_lab,
        show_auc = show_auc,
        tau = tau,
        title = title,
        x_breaks = x_breaks,
        x_lim = x_lim,
        x_name = x_name,
        y_breaks = y_breaks,
        y_lim = y_lim,
        y_name = y_name
      )
      return(out)
    }
    mcf_df <- .UniTwoArmMCFFrame(data = data_norm, tau = tau)
    out <- .UniTwoArmMCFPanel(
      mcf_df = mcf_df,
      arm_style = arm_style,
      show_auc = show_auc,
      tau = tau,
      title = title,
      x_breaks = x_breaks,
      x_lim = x_lim,
      x_name = x_name,
      y_breaks = y_breaks,
      y_lim = y_lim,
      y_name = y_name
    )
    return(out)
  }

  if (one_arm) {
    nar_df <- .UniOneArmNARFrame(
      data = data_norm,
      x_breaks = x_breaks,
      y_lab = color_lab
    )
    out <- .UniOneArmNARPanel(
      nar_df = nar_df,
      x_breaks = x_breaks,
      x_labs = x_labs,
      x_max = x_max,
      x_name = x_name
    )
    return(out)
  }

  nar_df <- .UniTwoArmNARFrame(
    data = data_norm,
    x_breaks = x_breaks,
    color_labs = arm_style$color_labs
  )
  out <- .UniTwoArmNARPanel(
    nar_df = nar_df,
    x_breaks = x_breaks,
    x_labs = x_labs,
    x_max = x_max,
    x_name = x_name
  )
  return(out)
}


# -----------------------------------------------------------------------------


#' Plot Area Under the Mean Cumulative Function
#'
#' @description
#' Deprecated: use \code{\link{PlotMCFs}} with \code{which_arm}, \code{show_auc =
#' TRUE}, and \code{show_nar = FALSE}.
#'
#' @param data Data.frame.
#' @param which_arm Arm to plot.
#' @param arm_label Label for the arm.
#' @param arm_name Name of arm column in data.
#' @param color Color.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param strata_name Name of stratum column in data.
#' @param tau Truncation time for shading.
#' @param time_name Name of time column in data.
#' @param title Plot title.
#' @param jump_weights Optional per-event jump weights.
#' @param x_breaks X-axis breaks.
#' @param x_lim X-axis limits.
#' @param x_name X-axis label.
#' @param y_breaks Y-axis breaks.
#' @param y_lim Y-axis limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' @importFrom dplyr "%>%"
#' @export
PlotAUMCF <- function(
    data,
    which_arm,
    arm_label = "Placebo",
    arm_name = "arm",
    color = "#C65842",
    idx_name = "idx",
    status_name = "status",
    strata_name = NULL,
    time_name = "time",
    title = NULL,
    tau = NULL,
    jump_weights = NULL,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  out <- PlotMCFs(
    data = data,
    which_arm = which_arm,
    arm_name = arm_name,
    color = color,
    color_lab = arm_label,
    idx_name = idx_name,
    status_name = status_name,
    strata_name = strata_name,
    time_name = time_name,
    title = title,
    tau = tau,
    jump_weights = jump_weights,
    x_breaks = x_breaks,
    x_lim = x_lim,
    x_name = x_name,
    y_breaks = y_breaks,
    y_lim = y_lim,
    y_name = y_name,
    show_nar = FALSE,
    show_auc = TRUE
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Two Sample Number at Risk Plotting Frame
#'
#' Two sample numbers at risk for recurrent events data.
#'
#' @param data Data.frame.
#' @param x_breaks Time points at which to determine the NARs.
#' @param arm_name Name of arm column.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame containing `time`, `nar_ctrl`, `nar_trt`.
#' @importFrom dplyr "%>%"
TwoSampleNARFrame <- function(
    data,
    x_breaks,
    arm_name = "arm",
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
) {

  df <- .NormDataForPlot(
    data = data,
    arm_name = arm_name,
    strata_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    jump_weights = NULL
  )

  arm <- NULL
  g0 <- df %>%
    dplyr::filter(arm == 0) %>%
    NARCurve()
  g1 <- df %>%
    dplyr::filter(arm == 1) %>%
    NARCurve()

  out <- data.frame(
    time = x_breaks,
    nar_ctrl = g0(x_breaks),
    nar_trt = g1(x_breaks)
  )
  return(out)
}


#' Plot Two Sample Number at Risk
#'
#' @description
#' Deprecated: use \code{\link{PlotMCFs}} with \code{show_mcf = FALSE} and
#' \code{show_nar = TRUE}.
#'
#' @param data Data.frame.
#' @param x_breaks X-axis breaks.
#' @param arm_name Name of arm column.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param x_labs X-axis tick labels.
#' @param x_name X-axis label.
#' @param x_max X-axis upper limit.
#' @param y_labs Y-axis tick labels.
#' @return ggplot.
#' @export
PlotNARs <- function(
    data,
    x_breaks,
    arm_name = "arm",
    idx_name = "idx",
    status_name = "status",
    time_name = "time",
    x_labs = NULL,
    x_max = NULL,
    x_name = NULL,
    y_labs = c("Ctrl", "Trt")
) {

  out <- PlotMCFs(
    data = data,
    arm_name = arm_name,
    color_labs = y_labs,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    x_breaks = x_breaks,
    x_labs = x_labs,
    x_max = x_max,
    x_name = x_name,
    show_mcf = FALSE,
    show_nar = TRUE
  )
  return(out)
}
