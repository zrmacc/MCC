# Purpose: Function to plot the mean cumulative functions,
# comparing two treatment arms.
# Updated: 2021-08-06

# -----------------------------------------------------------------------------

#' MCF Curve
#' 
#' Construct a function to evaluate the mean cumulative function at a given time
#' for a single treatment arm.
#'   
#' @param data Data.frame.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param time_name Name of column column in data.
#' @return stepfun.
#' 
#' @export

MCFCurve <- function(
  data, 
  idx_name = "idx",
  status_name = "status",
  time_name = "time"
) {
  
  # Data preparation.
  df <- data %>%
    dplyr::rename(
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
    )
  
  # Construct MCF.
  mcf <- MCC::CalcMCF(
    time = df$time, 
    status = df$status, 
    idx = df$idx,
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
#' 
#' @importFrom dplyr "%>%"
#' @export

NARCurve <- function(
  data, 
  idx_name = "idx",
  status_name = "status",
  time_name = "time"
) {
  
  # Data preparation.
  df <- data %>%
    dplyr::rename(
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
    )
  
  # Fit cumulative incidence curve.
  fit <- MCC::CalcMCF(data$time, data$status, data$idx, calc_var = FALSE)
  g <- stats::stepfun(x = fit$time, y = c(length(unique(df$idx)), fit$nar))
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
#' @return Data.frame containing `time` and `mcf`.

MCFPlotFrame <- function(
  data,
  eval_points = 1000,
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time"
) {
  
  # Data preparation.
  df <- data %>%
    dplyr::rename(
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
    )
  g <- df %>% MCFCurve()
  
  # Time grid.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(
    time = times,
    mcf = g(times)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Plot MCFs.
#' 
#' Plot the mean cumulative functions comparing two treatment arms.
#'
#' @param data Data.frame.
#' @param arm_name Name of arm column in data.
#' @param color_labs Color labels.
#' @param ctrl_color Color for control arm.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column in data.
#' @param tau Truncation time.
#' @param time_name Name of column column in data.
#' @param title Plot title.
#' @param trt_color Color for treatment arm.
#' @param x_breaks X-axis breaks.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis label.
#' @param y_lim Y-axis limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' 
#' @importFrom dplyr "%>%"
#' @importFrom ggplot2 aes
#' @export

PlotMCFs <- function(
  data,
  arm_name = "arm",
  color_labs = c("Ctrl", "Trt"),
  ctrl_color = "#C65842",
  idx_name = "idx",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  trt_color = "#6385B8",
  x_breaks = NULL,
  x_max = NULL,
  x_name = "Time",
  y_lim = NULL,
  y_name = "Mean Cumulative Count"
) {
  
  # Data preparation.
  data <- data %>%
    dplyr::rename(
      "arm" = {{arm_name}},
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
    )
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data$time)
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  
  # Split data.
  data0 <- data %>% dplyr::filter(arm == 0)
  data1 <- data %>% dplyr::filter(arm == 1)
  
  # Estimate mean cumulative function (MCF).
  fit_mcf_0 <- MCC::CalcMCF(
    time = data0$time,
    status = data0$status,
    idx = data0$idx
  )
  
  fit_mcf_1 <- MCC::CalcMCF(
    time = data1$time,
    status = data1$status,
    idx = data1$idx
  )
  
  # MCF function for arm 0
  g0 <- stats::stepfun(
    x = fit_mcf_0$time,
    y = c(0, fit_mcf_0$mcf)
  )
  
  # MCF function for arm 1
  g1 <- stats::stepfun(
    x = fit_mcf_1$time,
    y = c(0, fit_mcf_1$mcf)
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
      legend.position = c(0.2, 0.8)
    ) + 
    ggplot2::geom_step(
      data = df, 
      aes(x = time, y = mcf, color = arm), 
      size = 1) + 
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(ctrl_color, trt_color),
      labels = color_labs
    ) + 
    ggplot2::scale_x_continuous(
      name = x_name
    ) +
    ggplot2::scale_y_continuous(
      name = y_name,
      limits = y_lim
    ) + 
    ggplot2::ggtitle(
      label = title
    )
  
  # Output.
  return(q)
}


# -----------------------------------------------------------------------------

#' Plot AUMCFs.
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
#' @param tau Truncation time for shading.
#' @param time_name Name of column column in data.
#' @param title Plot title.
#' @param x_breaks X-axis breaks.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis label.
#' @param y_lim Y-axis limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' 
#' @importFrom dplyr "%>%"
#' @importFrom ggplot2 aes element_blank
#' @export

PlotAUMCFs <- function(
  data,
  which_arm,
  arm_label = "Placebo",
  arm_name = "arm",
  color = "#C65842",
  idx_name = "idx",
  status_name = "status",
  time_name = "time",
  title = NULL,
  tau = NULL,
  x_breaks = NULL,
  x_max = NULL,
  x_name = "Time",
  y_lim = NULL,
  y_name = "Mean Cumulative Count"
) {
  
  # Data preparation.
  data <- data %>%
    dplyr::rename(
      "arm" = {{arm_name}},
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
    )
  
  if (is.null(x_max)) {
    x_max <- max(data$time)
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  
  # Split data.
  arm <- NULL
  data <- data %>% dplyr::filter(arm == which_arm)
  
  # Estimate mean cumulative function (MCF).
  fit_mcf <- MCC::CalcMCF(
    time = data$time,
    status = data$status,
    idx = data$idx
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
  mcf <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.2, 0.8)
    ) + 
    ggplot2::geom_ribbon(
      data = df_shade,
      aes(x = time, ymin = 0, ymax = mcf, fill = arm),
      alpha = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = color,
      labels = arm_label
    ) +
    ggplot2::geom_step(
      data = df, 
      aes(x = time, y = mcf), 
      color = color,
      size = 1
    ) + 
    ggplot2::scale_x_continuous(
      name = x_name
    ) +
    ggplot2::scale_y_continuous(
      name = y_name,
      limits = y_lim
    ) + 
    ggplot2::ggtitle(
      label = title
    )
  
  # Output.
  return(q)
}


# -----------------------------------------------------------------------------

#' Number at Risk Plotting Frame
#' 
#' Numbers at risk for recurrent events data.
#' 
#' @param data Data.frame.
#' @param x_breaks Time points at which to determine the NARs.
#' @param arm_name Name of arm column.
#' @param idx_name Name of index (subject identifier) column in data.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame containing `time`, `nar_ctrl`, `nar_trt`.
#' 
#' @importFrom dplyr "%>%" filter rename 

NARPlotFrame <- function(
  data, 
  x_breaks, 
  arm_name = "arm",
  idx_name = "idx",
  status_name = "status",
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      "arm" = {{arm_name}},
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
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


#' Plot Numbers at Risk
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
#' 
#' @importFrom dplyr "%>%"
#' @importFrom ggplot2 aes element_blank
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
  
  # Data prep.
  nar_ctrl <- NULL
  nar_trt <- NULL
  df <- data %>%
    dplyr::rename(
      "arm" = {{arm_name}},
      "idx" = {{idx_name}},
      "status" = {{status_name}},
      "time" = {{time_name}}
    ) %>%
    NARPlotFrame(x_breaks = x_breaks) %>%
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
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ggplot2::geom_text(
      aes(x = time, y = arm, label = nar)
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

