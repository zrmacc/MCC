# Purpose: Function to plot the mean cumulative functions,
# comparing two treatment arms.
# Updated: 2026-02-04

# ------------------------------------------------------------------------------

#' Plot Two Sample Mean Cumulative Function
#' 
#' Plot the mean cumulative functions comparing two treatment arms.
#'
#' @param data Data.frame.
#' @param arm_name Name of arm column in data.
#' @param color_labs Color labels.
#' @param ctrl_color Color for control arm.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param strata_name Name of stratum column in data. 
#' @param tau Truncation time.
#' @param time_name Name of column column in data.
#' @param title Plot title.
#' @param trt_color Color for treatment arm.
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @param x_breaks X-axis breaks.
#' @param x_lim X-axis limits.
#' @param x_name X-axis label.
#' @param y_breaks Y-axis breaks.
#' @param y_lim Y-axis limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' @export
PlotMCFs <- function(
  data,
  arm_name = "arm",
  color_labs = c("Ctrl", "Trt"),
  ctrl_color = "#C65842",
  idx_name = "idx",
  status_name = "status",
  strata_name = NULL,
  tau = NULL,
  time_name = "time",
  title = NULL,
  trt_color = "#6385B8",
  weights = NULL,
  x_breaks = NULL,
  x_lim = NULL,
  x_name = "Time",
  y_breaks = NULL,
  y_lim = NULL,
  y_name = "Mean Cumulative Count"
) {
  
  # Data preparation.
  data <- .NormDataForPlot(
    data = data,
    arm_name = arm_name,
    strata_name = strata_name,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = weights
  )
  
  # Truncation.
  if (is.null(x_lim[2])) {
    x_max <- max(data$time)
  } else{
    x_max <- x_lim[2]
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  
  # Calculate marginal MCF.
  marg_mcf <- CalcMargMCF(data)
  
  # MCF function for arm 0
  g0 <- stats::stepfun(
    x = marg_mcf$time[marg_mcf$arm == 0],
    y = c(0, marg_mcf$mcf[marg_mcf$arm == 0])
  )
  
  # MCF function for arm 1
  g1 <- stats::stepfun(
    x = marg_mcf$time[marg_mcf$arm == 1],
    y = c(0, marg_mcf$mcf[marg_mcf$arm == 1])
  )
  
  # Plotting frame for control arm.
  df0 <- data.frame(time = seq(from = 0, to = x_max, length.out = 200))
  df1 <- df0
  df0$mcf <- g0(df0$time)
  df0$arm <- 0
  
  # Plotting frame for treatment arm.
  df1$mcf <- g1(df1$time)
  df1$arm <- 1
  
  df <- rbind(df0, df1)
  df$arm <- factor(df$arm, levels = c(0, 1))
  
  # Plotting.
  arm <- NULL
  mcf <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position.inside = c(0.2, 0.8)
    ) + 
    ggplot2::geom_step(
      data = df, 
      ggplot2::aes(x = time, y = mcf, color = arm), 
      linewidth = 1
    ) + 
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(ctrl_color, trt_color),
      labels = color_labs
    )
  
  # X-axis.
  if (is.null(x_breaks)) {
    q <- q + 
      ggplot2::scale_x_continuous(
        name = x_name,
        limits = x_lim
      )
  } else {
    q <- q + 
      ggplot2::scale_x_continuous(
        name = x_name,
        breaks = x_breaks,
        limits = x_lim
      )
  }

  # Y-axis.
  if (is.null(y_breaks)) {
    q <- q + 
      ggplot2::scale_y_continuous(
        name = y_name,
        limits = y_lim
      )
  } else {
    q <- q + 
      ggplot2::scale_y_continuous(
        name = y_name,
        breaks = y_breaks,
        limits = y_lim
      )
  }
  
  # Title.
  q <- q + 
    ggplot2::ggtitle(
      label = title
    )
  
  # Output.
  return(q)
}


# -----------------------------------------------------------------------------


#' Plot Area Under the Mean Cumulative Function
#' 
#' Plot area under the mean cumulative function for a single treatment arm.
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
#' @param time_name Name of column column in data.
#' @param title Plot title.
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
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
  weights = NULL,
  x_breaks = NULL,
  x_lim = NULL,
  x_name = "Time",
  y_breaks = NULL,
  y_lim = NULL,
  y_name = "Mean Cumulative Count"
) {
  
  # Data preparation.
  data <- .NormDataForPlot(
    data = data,
    arm_name = arm_name,
    strata_name = strata_name,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = weights
  )
  
  # Truncation.
  if (is.null(x_lim[2])) {
    x_max <- max(data$time)
  } else {
    x_max <- x_lim[2]
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  
  # Split data.
  arm <- NULL
  fit_mcf <- CalcMargMCF(data) %>% dplyr::filter(arm == which_arm)
  
  # MCF function.
  g <- stats::stepfun(
    x = fit_mcf$time,
    y = c(0, fit_mcf$mcf)
  )
  
  # Plotting frames.
  df <- data.frame(time = seq(from = 0, to = x_max, length.out = 1001))
  df$mcf <- g(df$time)
  df$arm <- 0
  df_shade <- df %>% dplyr::filter(time <= tau)
  df_shade$arm <- factor(df_shade$arm)
  
  # Plotting.
  mcf <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position.inside = c(0.2, 0.8)
    ) + 
    ggplot2::geom_ribbon(
      data = df_shade,
      ggplot2::aes(x = time, ymin = 0, ymax = mcf, fill = arm),
      alpha = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = color,
      labels = arm_label
    ) +
    ggplot2::geom_step(
      data = df, 
      ggplot2::aes(x = time, y = mcf), 
      color = color,
      linewidth = 1
    ) 
  
  # X-axis.
  if (is.null(x_breaks)) {
      q <- q + 
        ggplot2::scale_x_continuous(
          name = x_name,
          limits = x_lim
        )
    } else {
      q <- q + 
        ggplot2::scale_x_continuous(
          name = x_name,
          breaks = x_breaks,
          limits = x_lim
        )
    }
  
  # Y-axis.
  if (is.null(y_breaks)) {
    q <- q + 
      ggplot2::scale_y_continuous(
        name = y_name,
        limits = y_lim
      )
  } else {
    q <- q + 
      ggplot2::scale_y_continuous(
        name = y_name,
        breaks = y_breaks,
        limits = y_lim
      )
  }
  
  # Title.
  q <- q + 
    ggplot2::ggtitle(
      label = title
    )
  
  # Output.
  return(q)
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
  
  # Data preparation.
  df <- .NormDataForPlot(
    data = data,
    arm_name = arm_name,
    strata_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = NULL
  )
  
  # NAR functions.
  arm <- NULL
  g0 <- df %>% dplyr::filter(arm == 0) %>% MCC::NARCurve()
  g1 <- df %>% dplyr::filter(arm == 1) %>% MCC::NARCurve()
  
  # Output.
  out <- data.frame(
    time = x_breaks,
    nar_ctrl = g0(x_breaks),
    nar_trt = g1(x_breaks)
  )
  return(out)
}


#' Plot Two Sample Number at Risk
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
  
  # Defaults.
  if (is.null(x_labs)) {
    x_labs = x_breaks
  }
  if (is.null(x_max)) {
    x_max = max(x_breaks)
  }
  
  # Data preparation.
  df <- .NormDataForPlot(
    data = data,
    arm_name = arm_name,
    strata_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = NULL
  )
  
  nar_ctrl <- NULL
  nar_trt <- NULL
  df <- TwoSampleNARFrame(df, x_breaks = x_breaks) %>%
    tidyr::pivot_longer(
      cols = c(nar_ctrl, nar_trt),
      names_to = "arm",
      values_to = "nar"
    ) %>%
    dplyr::mutate(
      arm = factor(arm, c("nar_ctrl", "nar_trt"), y_labs)
    )
  
  # Plotting.
  arm <- NULL
  nar <- NULL
  time <- NULL
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = time, y = arm, label = nar)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) + 
    ggplot2::scale_y_discrete(
      name = NULL,
      labels = y_labs
    )
  return(q)
}
