# Purpose: Function to plot the mean cumulative function for one treatment arm.
# Updated: 2026-06-17

# ------------------------------------------------------------------------------


#' MCF Curve
#' 
#' Construct a function to evaluate the mean cumulative function at a given time
#' for a single treatment arm.
#'   
#' @param data Data.frame.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param time_name Name of column column in data.
#' @param jump_weights Optional column of jump weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return stepfun.
#' @export
MCFCurve <- function(
  data, 
  idx_name = "idx",
  status_name = "status",
  time_name = "time",
  jump_weights = NULL
) {
  
  df <- .NormDataForPlot(
    data = data,
    arm_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    jump_weights = jump_weights
  )
  mcf <- CalcMCF(
    idx = df$idx,
    status = df$status, 
    time = df$time, 
    jump_weights = df$jump_weights,
    calc_var = FALSE
  )
  
  g <- stats::stepfun(x = mcf$time, y = c(0, mcf$mcf))
  return(g)
}


# -----------------------------------------------------------------------------


#' Number at Risk Curve
#' 
#' Return a function that calculates the number at risk for a single
#' treatment arm.
#' 
#' @param data Data.frame.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param time_name Name of column column in data.
#' @return stepfun.
#' @export
NARCurve <- function(
  data, 
  idx_name = "idx",
  status_name = "status",
  time_name = "time"
) {
  
  # Data preparation.
  df <- .NormDataForPlot(
    data = data,
    arm_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    jump_weights = NULL
  )
  
  # Fit cumulative incidence curve.
  fit <- CalcMCF(
    idx = df$idx,
    status = df$status,
    time = df$time,
    jump_weights = df$jump_weights,
    calc_var = FALSE
  )
  
  # Case where last observation is censoring or death.
  last_row <- fit[nrow(fit), ]
  last_row$time <- last_row$time + 1e-8
  last_row$nar <- last_row$nar - (last_row$censor + last_row$death)
  fit <- rbind(fit, last_row)
  
  g <- stats::stepfun(
    x = fit$time, 
    y = c(length(unique(df$idx)), fit$nar)
  )
  return(g)
}


# -----------------------------------------------------------------------------

#' MCF Plotting Frame
#' 
#' Constructs the MCF plotting frame for a single treatment arm.
#' 
#' @param data Data.frame.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param tau Truncation time.
#' @param time_name Name of column column in data.
#' @param jump_weights Optional column of jump weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Data.frame containing `time` and `mcf`.
MCFPlotFrame <- function(
  data,
  eval_points = 1000,
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  jump_weights = NULL
) {
  
  # Data preparation.
  df <- .NormDataForPlot(
    data = data,
    arm_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    jump_weights = jump_weights
  )
  
  # MCF curve.
  g <- MCFCurve(data = df, jump_weights = df$jump_weights)
  
  # Time grid.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(time = times, mcf = g(times))
  return(out)
}


# -----------------------------------------------------------------------------
# Internal univariate plotting helpers.
# -----------------------------------------------------------------------------

#' ggplot theme for univariate MCF panels.
#'
#' @return ggplot2 theme object.
#' @noRd
.UniGgTheme <- function() {

  out <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
  return(out)
}


#' ggplot theme for univariate NAR panels.
#'
#' @return ggplot2 theme object.
#' @noRd
.UniNARTheme <- function() {

  out <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  return(out)
}


#' Resolve two-arm colors and labels for univariate plots.
#'
#' @param color_labs Arm labels.
#' @param colors Arm colors.
#' @param ctrl_color Deprecated control color.
#' @param trt_color Deprecated treatment color.
#' @return List with color_labs and colors vectors of length 2.
#' @noRd
.UniResolveArmStyle <- function(
    color_labs = c("Ctrl", "Trt"),
    colors = NULL,
    ctrl_color = "#C65842",
    trt_color = "#6385B8"
) {

  lab_vec <- as.character(color_labs)
  if (length(lab_vec) != 2L) {
    stop("color_labs must have length 2 for two-arm plots.", call. = FALSE)
  }
  if (is.null(colors)) {
    col_vec <- c(ctrl_color, trt_color)
  } else {
    col_vec <- as.character(colors)
    if (length(col_vec) < 2L) {
      stop("colors must have at least two entries for two-arm plots.", call. = FALSE)
    }
    col_vec <- col_vec[1:2]
  }
  out <- list(
    color_labs = lab_vec,
    colors = col_vec
  )
  return(out)
}


#' Infer one-arm vs two-arm layout from plotting data.
#'
#' @param data Data.frame after optional arm subsetting.
#' @param arm_name Name of arm column.
#' @return Logical scalar.
#' @noRd
.UniIsOneArm <- function(data, arm_name = "arm") {

  if (!arm_name %in% names(data)) {
    out <- TRUE
    return(out)
  }
  arms <- unique(data[[arm_name]])
  out <- length(arms) <= 1L
  return(out)
}


#' Apply continuous x/y scales to a univariate ggplot.
#'
#' @param q ggplot object.
#' @param x_breaks X breaks.
#' @param x_lim X limits.
#' @param x_name X-axis label.
#' @param y_breaks Y breaks.
#' @param y_lim Y limits.
#' @param y_name Y-axis label.
#' @return Updated ggplot object.
#' @noRd
.UniApplyAxes <- function(
    q,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  if (is.null(x_breaks)) {
    q <- q + ggplot2::scale_x_continuous(name = x_name, limits = x_lim)
  } else {
    q <- q + ggplot2::scale_x_continuous(
      name = x_name,
      breaks = x_breaks,
      limits = x_lim
    )
  }
  if (is.null(y_breaks)) {
    q <- q + ggplot2::scale_y_continuous(name = y_name, limits = y_lim)
  } else {
    q <- q + ggplot2::scale_y_continuous(
      name = y_name,
      breaks = y_breaks,
      limits = y_lim
    )
  }
  return(q)
}


#' One-arm MCF plotting frame.
#'
#' @param data Formatted one-arm data.frame.
#' @param tau Truncation time.
#' @param eval_points Number of evaluation points.
#' @return Data.frame with time and mcf.
#' @noRd
.UniOneArmMCFFrame <- function(data, tau, eval_points = 200L) {

  mcf <- CalcMCF(
    idx = data[["idx"]],
    status = data[["status"]],
    time = data[["time"]],
    jump_weights = data[["jump_weights"]],
    calc_var = FALSE
  )
  g <- stats::stepfun(x = mcf[["time"]], y = c(0, mcf[["mcf"]]))
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(time = times, mcf = g(times))
  return(out)
}


#' Two-arm MCF plotting frame.
#'
#' @param data Formatted two-arm data.frame.
#' @param tau Truncation time.
#' @param eval_points Number of evaluation points per arm.
#' @return Data.frame with time, mcf, and arm factor.
#' @noRd
.UniTwoArmMCFFrame <- function(data, tau, eval_points = 200L) {

  marg_mcf <- CalcMargMCF(data)
  g0 <- stats::stepfun(
    x = marg_mcf[["time"]][marg_mcf[["arm"]] == 0],
    y = c(0, marg_mcf[["mcf"]][marg_mcf[["arm"]] == 0])
  )
  g1 <- stats::stepfun(
    x = marg_mcf[["time"]][marg_mcf[["arm"]] == 1],
    y = c(0, marg_mcf[["mcf"]][marg_mcf[["arm"]] == 1])
  )
  times <- seq(from = 0, to = tau, length.out = eval_points)
  df0 <- data.frame(time = times, mcf = g0(times), arm = 0)
  df1 <- data.frame(time = times, mcf = g1(times), arm = 1)
  out <- rbind(df0, df1)
  rownames(out) <- NULL
  out[["arm"]] <- factor(out[["arm"]], levels = c(0, 1))
  return(out)
}


#' One-arm NAR plotting frame.
#'
#' @param data Formatted one-arm data.frame.
#' @param x_breaks Time points for NAR evaluation.
#' @param y_lab Y-axis row label.
#' @return Data.frame with time, nar, and arm factor.
#' @noRd
.UniOneArmNARFrame <- function(data, x_breaks, y_lab = "Arm") {

  g <- NARCurve(data = data)
  out <- data.frame(
    time = x_breaks,
    nar = g(x_breaks)
  )
  out[["arm"]] <- factor(y_lab, levels = y_lab)
  return(out)
}


#' Two-arm NAR plotting frame.
#'
#' @param data Formatted two-arm data.frame.
#' @param x_breaks Time points for NAR evaluation.
#' @param color_labs Arm labels for the y-axis.
#' @return Data.frame with time, nar, and arm factor.
#' @noRd
.UniTwoArmNARFrame <- function(data, x_breaks, color_labs = c("Ctrl", "Trt")) {

  arm <- NULL
  g0 <- data %>%
    dplyr::filter(arm == 0) %>%
    NARCurve()
  g1 <- data %>%
    dplyr::filter(arm == 1) %>%
    NARCurve()
  out <- data.frame(
    time = rep(x_breaks, 2L),
    nar = c(g0(x_breaks), g1(x_breaks)),
    arm = factor(
      rep(c(0L, 1L), each = length(x_breaks)),
      levels = c(0, 1),
      labels = color_labs
    )
  )
  return(out)
}


#' One-arm MCF panel.
#'
#' @param mcf_df MCF data.frame.
#' @param color Line and fill color.
#' @param color_lab Legend label.
#' @param show_auc Shade area under the curve?
#' @param tau Truncation time for shading.
#' @param title Panel title.
#' @param x_breaks X breaks.
#' @param x_lim X limits.
#' @param x_name X-axis label.
#' @param y_breaks Y breaks.
#' @param y_lim Y limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' @noRd
.UniOneArmMCFPanel <- function(
    mcf_df,
    color = "#C65842",
    color_lab = "Arm",
    show_auc = FALSE,
    tau = NULL,
    title = NULL,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  mcf <- time <- NULL
  q <- ggplot2::ggplot() +
    .UniGgTheme()
  if (show_auc) {
    mcf_shade <- mcf_df[mcf_df[["time"]] <= tau, , drop = FALSE]
    q <- q +
      ggplot2::geom_ribbon(
        data = mcf_shade,
        ggplot2::aes(x = time, ymin = 0, ymax = mcf),
        fill = color,
        alpha = 0.5
      )
  }
  q <- q +
    ggplot2::geom_step(
      data = mcf_df,
      ggplot2::aes(x = time, y = mcf, color = "curve"),
      linewidth = 1
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(curve = color),
      labels = color_lab
    )
  q <- .UniApplyAxes(
    q = q,
    x_breaks = x_breaks,
    x_lim = x_lim,
    x_name = x_name,
    y_breaks = y_breaks,
    y_lim = y_lim,
    y_name = y_name
  )
  q <- q + ggplot2::ggtitle(label = title)
  return(q)
}


#' Two-arm MCF panel.
#'
#' @param mcf_df MCF data.frame with arm factor.
#' @param arm_style List from `.UniResolveArmStyle`.
#' @param show_auc Shade area under curves?
#' @param tau Truncation time for shading.
#' @param title Panel title.
#' @param x_breaks X breaks.
#' @param x_lim X limits.
#' @param x_name X-axis label.
#' @param y_breaks Y breaks.
#' @param y_lim Y limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' @noRd
.UniTwoArmMCFPanel <- function(
    mcf_df,
    arm_style,
    show_auc = FALSE,
    tau = NULL,
    title = NULL,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  arm <- mcf <- time <- NULL
  q <- ggplot2::ggplot() +
    .UniGgTheme()
  if (show_auc) {
    mcf_shade <- mcf_df[mcf_df[["time"]] <= tau, , drop = FALSE]
    q <- q +
      ggplot2::geom_ribbon(
        data = mcf_shade,
        ggplot2::aes(x = time, ymin = 0, ymax = mcf, fill = arm),
        alpha = 0.5
      ) +
      ggplot2::scale_fill_manual(
        values = arm_style$colors,
        guide = "none"
      )
  }
  q <- q +
    ggplot2::geom_step(
      data = mcf_df,
      ggplot2::aes(x = time, y = mcf, color = arm),
      linewidth = 1
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = arm_style$colors,
      labels = arm_style$color_labs
    )
  q <- .UniApplyAxes(
    q = q,
    x_breaks = x_breaks,
    x_lim = x_lim,
    x_name = x_name,
    y_breaks = y_breaks,
    y_lim = y_lim,
    y_name = y_name
  )
  q <- q + ggplot2::ggtitle(label = title)
  return(q)
}


#' One-arm NAR panel.
#'
#' @param nar_df NAR data.frame.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis tick labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis label.
#' @return ggplot object.
#' @noRd
.UniOneArmNARPanel <- function(
    nar_df,
    x_breaks,
    x_labs = NULL,
    x_max = NULL,
    x_name = NULL
) {

  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  if (is.null(x_max)) {
    x_max <- max(x_breaks)
  }
  arm <- nar <- time <- NULL
  q <- ggplot2::ggplot(data = nar_df) +
    .UniNARTheme() +
    ggplot2::geom_text(ggplot2::aes(x = time, y = arm, label = nar)) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_discrete(name = NULL)
  return(q)
}


#' Two-arm NAR panel.
#'
#' @param nar_df NAR data.frame.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis tick labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis label.
#' @return ggplot object.
#' @noRd
.UniTwoArmNARPanel <- function(
    nar_df,
    x_breaks,
    x_labs = NULL,
    x_max = NULL,
    x_name = NULL
) {

  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  if (is.null(x_max)) {
    x_max <- max(x_breaks)
  }
  arm <- nar <- time <- NULL
  q <- ggplot2::ggplot(data = nar_df) +
    .UniNARTheme() +
    ggplot2::geom_text(ggplot2::aes(x = time, y = arm, label = nar)) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_discrete(name = NULL)
  return(q)
}


#' Stack univariate MCF and NAR panels.
#'
#' @param one_arm Logical; one-arm vs two-arm layout.
#' @param data_norm Formatted plotting data.frame.
#' @param arm_style Two-arm style list or NULL.
#' @param tau Truncation time.
#' @param x_breaks NAR evaluation times.
#' @param x_labs X-axis tick labels for NAR row.
#' @param x_max X-axis upper limit for NAR row.
#' @param color Line color for one-arm plots.
#' @param color_lab Legend label for one-arm plots.
#' @param show_auc Shade area under MCF curves?
#' @param title Plot title.
#' @param x_lim X limits.
#' @param x_name X-axis label for NAR row.
#' @param y_breaks Y breaks for MCF row.
#' @param y_lim Y limits for MCF row.
#' @param y_name Y-axis label for MCF row.
#' @return cowplot grid object.
#' @noRd
.UniBuildCombinedPlot <- function(
    one_arm,
    data_norm,
    arm_style,
    tau,
    x_breaks,
    x_labs = NULL,
    x_max = NULL,
    color = "#C65842",
    color_lab = "Arm",
    show_auc = FALSE,
    title = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  if (is.null(x_max)) {
    x_max <- max(x_breaks)
  }
  if (!is.null(x_lim)) {
    x_max <- x_lim[2L]
  }

  if (one_arm) {
    mcf_df <- .UniOneArmMCFFrame(data = data_norm, tau = tau)
    nar_df <- .UniOneArmNARFrame(
      data = data_norm,
      x_breaks = x_breaks,
      y_lab = color_lab
    )
    p_mcf <- .UniOneArmMCFPanel(
      mcf_df = mcf_df,
      color = color,
      color_lab = color_lab,
      show_auc = show_auc,
      tau = tau,
      x_breaks = x_breaks,
      x_lim = x_lim,
      x_name = NULL,
      y_breaks = y_breaks,
      y_lim = y_lim,
      y_name = y_name
    )
    p_nar <- .UniOneArmNARPanel(
      nar_df = nar_df,
      x_breaks = x_breaks,
      x_labs = x_labs,
      x_max = x_max,
      x_name = x_name
    )
  } else {
    mcf_df <- .UniTwoArmMCFFrame(data = data_norm, tau = tau)
    nar_df <- .UniTwoArmNARFrame(
      data = data_norm,
      x_breaks = x_breaks,
      color_labs = arm_style$color_labs
    )
    p_mcf <- .UniTwoArmMCFPanel(
      mcf_df = mcf_df,
      arm_style = arm_style,
      show_auc = show_auc,
      tau = tau,
      x_breaks = x_breaks,
      x_lim = x_lim,
      x_name = NULL,
      y_breaks = y_breaks,
      y_lim = y_lim,
      y_name = y_name
    )
    p_nar <- .UniTwoArmNARPanel(
      nar_df = nar_df,
      x_breaks = x_breaks,
      x_labs = x_labs,
      x_max = x_max,
      x_name = x_name
    )
  }

  out <- cowplot::plot_grid(
    p_mcf,
    p_nar,
    ncol = 1L,
    rel_heights = c(3, 1),
    align = "v",
    axis = "l"
  )
  if (!is.null(title)) {
    out <- cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label(title, fontface = "bold"),
      out,
      ncol = 1,
      rel_heights = c(0.05, 0.95)
    )
  }
  return(out)
}


# -----------------------------------------------------------------------------

#' Plot One Sample Mean Cumulative Function
#'
#' @description
#' Deprecated: use \code{\link{PlotMCFs}} with \code{show_nar = FALSE} and
#' \code{show_auc = FALSE}.
#'
#' @inheritParams PlotMCFs
#' @return ggplot object.
#' @importFrom dplyr "%>%"
#' @export
PlotOneSampleMCF <- function(
  data,
  color = "#C65842",
  color_lab = "Arm",
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
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
    color = color,
    color_lab = color_lab,
    idx_name = idx_name,
    status_name = status_name,
    tau = tau,
    time_name = time_name,
    title = title,
    jump_weights = jump_weights,
    x_breaks = x_breaks,
    x_lim = x_lim,
    x_name = x_name,
    y_breaks = y_breaks,
    y_lim = y_lim,
    y_name = y_name,
    show_nar = FALSE,
    show_auc = FALSE
  )
  return(out)
}


# -----------------------------------------------------------------------------


#' Plot One Sample Area Under the Mean Cumulative Function
#'
#' @description
#' Deprecated: use \code{\link{PlotMCFs}} with \code{show_nar = FALSE} and
#' \code{show_auc = TRUE}.
#'
#' @inheritParams PlotMCFs
#' @return ggplot object.
#' @importFrom dplyr "%>%"
#' @export
PlotOneSampleAUMCF <- function(
  data,
  color = "#C65842",
  color_lab = "Arm",
  idx_name = "idx",
  status_name = "status",
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
    color = color,
    color_lab = color_lab,
    idx_name = idx_name,
    status_name = status_name,
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

#' Plot One Sample Number at Risk
#'
#' @description
#' Deprecated: use \code{\link{PlotMCFs}} with \code{show_mcf = FALSE} and
#' \code{show_nar = TRUE}.
#'
#' @param data Data.frame.
#' @param x_breaks X-axis breaks.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param x_labs X-axis tick labels.
#' @param x_name X-axis label.
#' @param x_max X-axis upper limit.
#' @param y_lab Y-axis tick label.
#' @return ggplot.
#' @export
PlotOneSampleNAR <- function(
  data,
  x_breaks,
  idx_name = "idx",
  status_name = "status",
  time_name = "time",
  x_labs = NULL,
  x_max = NULL,
  x_name = NULL,
  y_lab = "Arm"
) {

  out <- PlotMCFs(
    data = data,
    color_lab = y_lab,
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

