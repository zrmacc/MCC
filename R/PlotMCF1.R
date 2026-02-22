# Purpose: Function to plot the mean cumulative function for one treatment arm.
# Updated: 2026-02-04

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
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return stepfun.
#' @export
MCFCurve <- function(
  data, 
  idx_name = "idx",
  status_name = "status",
  time_name = "time",
  weights = NULL
) {
  
  df <- .NormDataForPlot(
    data = data,
    arm_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = weights
  )
  mcf <- CalcMCF(
    idx = df$idx,
    status = df$status, 
    time = df$time, 
    weights = df$weights,
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
    weights = NULL
  )
  
  # Fit cumulative incidence curve.
  fit <- CalcMCF(
    idx = df$idx,
    status = df$status,
    time = df$time,
    weights = df$weights,
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
#' @param weights Optional column of weights, controlling the size of the jump
#'   in the cumulative count curve at times with status == 1.
#' @return Data.frame containing `time` and `mcf`.
MCFPlotFrame <- function(
  data,
  eval_points = 1000,
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  weights = NULL
) {
  
  # Data preparation.
  df <- .NormDataForPlot(
    data = data,
    arm_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = weights
  )
  
  # MCF curve.
  g <- MCFCurve(data = df, weights = df$weights)
  
  # Time grid.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(time = times, mcf = g(times))
  return(out)
}


# -----------------------------------------------------------------------------

#' Plot One Sample Mean Cumulative Function
#' 
#' Plot the mean cumulative function for a single treatment arm.
#'
#' @param data Data.frame.
#' @param color Color.
#' @param color_lab Color label.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param tau Truncation time.
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
PlotOneSampleMCF <- function(
  data,
  color = "#C65842",
  color_lab = "Arm",
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
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
    arm_name = NULL,
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
  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    weights = data$weights
  )
  
  # MCF function.
  g <- stats::stepfun(
    x = mcf$time,
    y = c(0, mcf$mcf)
  )
  
  # Plotting frame for control arm.
  df <- data.frame(time = seq(from = 0, to = x_max, length.out = 200))
  df$mcf <- g(df$time)
  
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
    ggplot2::geom_step(
      data = df, 
      ggplot2::aes(x = time, y = mcf, color = color), 
      linewidth = 1
    ) + 
    ggplot2::scale_color_manual(
      name = NULL,
      values = color,
      labels = color_lab
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


#' Plot One Sample Area Under the Mean Cumulative Function
#' 
#' Plot area under the mean cumulative function for a single treatment arm.
#'
#' @param data Data.frame.
#' @param color Color.
#' @param color_lab Color label.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
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
PlotOneSampleAUMCF <- function(
  data,
  color = "#C65842",
  color_lab = "Arm",
  idx_name = "idx",
  status_name = "status",
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
    arm_name = NULL,
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
  
  # Estimate mean cumulative function.
  fit_mcf <- MCC::CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    weights = data$weights
  )
  
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
  arm <- mcf <- time <- NULL
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
      labels = color_lab
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

#' Plot One Sample Number at Risk
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
    arm_name = NULL,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name,
    weights = NULL
  )
  
  # NAR data frame.
  g <- NARCurve(data = df)
  
  # Output.
  df <- data.frame(
    time = x_breaks,
    nar = g(x_breaks)
  )
  
  # Plotting.
  arm <- nar <- time <- NULL
  df$arm <- factor(x = 0, levels = 0, labels = y_lab)
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
      labels = y_lab
    )
  return(q)
}

