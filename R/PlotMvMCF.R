# Purpose: Plot mean cumulative functions for multivariate recurrent events.
# Updated: 2026-06-17

# -----------------------------------------------------------------------------
# Internal helpers.
# -----------------------------------------------------------------------------

#' Normalize multivariate data for plotting.
#'
#' @param data Data.frame in long format.
#' @param arm_name Name of arm column, or NULL for one-sample.
#' @param cens_after_last Censor after last recurrent event if no terminal row?
#' @param event_type_name Name of event type column.
#' @param idx_name Name of subject index column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return List with formatted data, event_types, K, univariate_path,
#'   cens_after_last.
#' @noRd
.NormMvDataForPlot <- function(
    data,
    arm_name = "arm",
    cens_after_last = TRUE,
    event_type_name = "event_type",
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
) {

  if (is.null(arm_name)) {
    plot_data <- data
    if (!"arm" %in% names(plot_data)) {
      plot_data$arm <- 0L
    }
    fmt <- FormatMvData(
      data = plot_data,
      arm_name = "arm",
      cens_after_last = cens_after_last,
      event_type_name = event_type_name,
      idx_name = idx_name,
      status_name = status_name,
      time_name = time_name
    )
    arms <- unique(fmt$data[["arm"]])
    if (length(arms) > 1L) {
      stop(
        "One-sample multivariate plots require a single arm. ",
        "Use which_arm or subset to one arm.",
        call. = FALSE
      )
    }
  } else {
    fmt <- FormatMvData(
      data = data,
      arm_name = arm_name,
      cens_after_last = cens_after_last,
      event_type_name = event_type_name,
      idx_name = idx_name,
      status_name = status_name,
      time_name = time_name
    )
  }

  out <- list(
    data = fmt$data,
    event_types = fmt$event_types,
    K = fmt$K,
    univariate_path = fmt$univariate_path,
    cens_after_last = if (fmt$univariate_path) {
      cens_after_last
    } else {
      fmt$cens_after_last
    }
  )
  return(out)
}


#' Subset and format one event process for plotting.
#'
#' @param data Formatted multivariate data.frame.
#' @param k Event type.
#' @param arm Arm to retain, or NULL for all arms.
#' @param cens_after_last Censor after last event if no terminal row?
#' @return Per-process data.frame with idx, status, time, jump_weights.
#' @noRd
.ProcessDataForPlotK <- function(
    data,
    k,
    arm = NULL,
    cens_after_last = TRUE
) {

  plot_data <- data
  if (!is.null(arm)) {
    plot_data <- plot_data[plot_data[["arm"]] == arm, , drop = FALSE]
  }

  at_risk_idx <- .MvAtRiskIdx(data = plot_data, k = k)
  plot_data <- plot_data[plot_data[["idx"]] %in% at_risk_idx, , drop = FALSE]

  proc <- .SubsetProcessK(
    data = plot_data,
    k = k,
    cens_after_last = cens_after_last
  )
  return(proc)
}


#' Resolve display labels for event types.
#'
#' @param event_types Integer vector of event type codes.
#' @param event_type_labs Optional character vector of labels, in the same order
#'   as \code{event_types}, or a named vector keyed by event type.
#' @return Named character vector mapping event type code to label.
#' @noRd
.MvResolveEventTypeLabs <- function(event_types, event_type_labs = NULL) {

  event_types <- as.integer(event_types)
  type_chr <- as.character(event_types)
  if (is.null(event_type_labs)) {
    out <- stats::setNames(paste0("k = ", event_types), type_chr)
    return(out)
  }
  if (!is.null(names(event_type_labs))) {
    names_chr <- names(event_type_labs)
    if (any(!names_chr %in% type_chr)) {
      stop(
        "event_type_labs names must match event types in the data.",
        call. = FALSE
      )
    }
    missing <- setdiff(type_chr, names_chr)
    if (length(missing) > 0L) {
      stop(
        "event_type_labs is missing label(s) for event type(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    out <- stats::setNames(
      as.character(event_type_labs[type_chr]),
      type_chr
    )
    return(out)
  }
  if (length(event_type_labs) != length(event_types)) {
    stop(
      "event_type_labs must have one label per event type.",
      call. = FALSE
    )
  }
  out <- stats::setNames(as.character(event_type_labs), type_chr)
  return(out)
}


#' Attach resolved event-type labels to multivariate plot prep.
#'
#' @param mv_prep Output of `.NormMvDataForPlot`.
#' @param event_type_labs Optional event-type labels.
#' @return Updated \code{mv_prep} list.
#' @noRd
.MvPrepEventTypeLabs <- function(mv_prep, event_type_labs = NULL) {

  mv_prep$event_type_labs <- .MvResolveEventTypeLabs(
    event_types = mv_prep$event_types,
    event_type_labs = event_type_labs
  )
  return(mv_prep)
}


#' Panel title for one event type.
#'
#' @param mv_prep Output of `.NormMvDataForPlot` with \code{event_type_labs}.
#' @param k Event type code.
#' @return Character scalar label.
#' @noRd
.MvEventTypePanelTitle <- function(mv_prep, k) {

  out <- mv_prep$event_type_labs[[as.character(k)]]
  return(out)
}


#' Unique sorted arm codes in plotting data.
#'
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @return Sorted vector of unique arm values.
#' @noRd
.MvUniqueArms <- function(data, arm_name = "arm") {

  if (!arm_name %in% names(data)) {
    out <- 0L
    return(out)
  }
  out <- sort(unique(data[[arm_name]]))
  return(out)
}


#' Default line colors for multi-arm plots.
#'
#' @param n Number of arms.
#' @return Character vector of hex colors.
#' @noRd
.MvDefaultArmColors <- function(n) {

  base <- c("#C65842", "#6385B8")
  if (n <= 2L) {
    out <- base[seq_len(n)]
    return(out)
  }
  extra <- grDevices::colorRampPalette(c("#5BA85B", "#9B6BB8", "#E8A838"))(n - 2L)
  out <- c(base, extra)
  return(out)
}


#' Resolve arm labels and colors.
#'
#' @param arms Sorted vector of arm codes.
#' @param color_labs Optional arm labels.
#' @param colors Optional arm colors.
#' @return List with arms, color_labs, and colors (named by arm).
#' @noRd
.MvResolveArmStyle <- function(arms, color_labs = NULL, colors = NULL) {

  n <- length(arms)
  arm_chr <- as.character(arms)
  if (is.null(color_labs)) {
    lab_vec <- paste0("Arm ", arms)
  } else {
    lab_vec <- as.character(color_labs)
    if (length(lab_vec) != n) {
      stop("color_labs must have one label per arm.", call. = FALSE)
    }
  }
  if (is.null(colors)) {
    col_vec <- .MvDefaultArmColors(n)
  } else {
    col_vec <- as.character(colors)
    if (length(col_vec) < n) {
      stop("colors must have at least one entry per arm.", call. = FALSE)
    }
    col_vec <- col_vec[seq_len(n)]
  }
  out <- list(
    arms = arms,
    color_labs = stats::setNames(lab_vec, arm_chr),
    colors = stats::setNames(col_vec, arm_chr)
  )
  return(out)
}


#' Subset data to a single arm.
#'
#' @param data Data.frame.
#' @param which_arm Arm code to retain, or NULL.
#' @param arm_name Name of arm column.
#' @return Subset data.frame.
#' @noRd
.MvSubsetWhichArm <- function(data, which_arm, arm_name = "arm") {

  if (is.null(which_arm)) {
    return(data)
  }
  if (!arm_name %in% names(data)) {
    stop("which_arm requires column ", arm_name, ".", call. = FALSE)
  }
  keep <- data[[arm_name]] == which_arm
  if (!any(keep)) {
    stop(
      "which_arm ",
      which_arm,
      " not found in ",
      arm_name,
      ".",
      call. = FALSE
    )
  }
  out <- data[keep, , drop = FALSE]
  return(out)
}


#' Add facet labels for event type.
#'
#' @param data Data.frame with event_type column.
#' @param event_type_labs Named character vector of labels, or \code{NULL} for
#'   the default \code{k = <type>} labels.
#' @return Data.frame with event_type_label factor.
#' @noRd
.MvAddEventTypeLabel <- function(data, event_type_labs = NULL) {

  et <- data[["event_type"]]
  types <- sort(unique(et))
  type_chr <- as.character(types)
  if (is.null(event_type_labs)) {
    label_vec <- paste0("k = ", types)
    row_labels <- paste0("k = ", et)
  } else {
    label_vec <- event_type_labs[type_chr]
    if (any(is.na(label_vec))) {
      stop("Missing event type label(s).", call. = FALSE)
    }
    row_labels <- event_type_labs[as.character(et)]
  }
  data[["event_type_label"]] <- factor(
    row_labels,
    levels = unname(label_vec)
  )
  return(data)
}


#' Shared ggplot theme for multivariate MCF plots.
#'
#' @return ggplot2 theme object.
#' @noRd
.MvGgTheme <- function() {

  out <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain")
    )
  return(out)
}


#' ggplot theme for multivariate NAR panels.
#'
#' @return ggplot2 theme object.
#' @noRd
.MvNARTheme <- function() {

  out <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain")
    )
  return(out)
}


#' Default x-axis breaks for multivariate NAR tables.
#'
#' @param tau Truncation or upper time limit.
#' @return Numeric vector of breaks.
#' @noRd
.MvDefaultXBreaks <- function(tau) {

  if (is.null(tau)) {
    stop("tau is required to derive x_breaks.", call. = FALSE)
  }
  if (abs(tau - round(tau)) < 1e-8) {
    out <- seq(0, as.integer(round(tau)), by = 1L)
  } else {
    out <- pretty(c(0, tau), n = 5L)
  }
  return(out)
}


#' One-sample multivariate NAR plotting frame.
#'
#' @param mv_prep Output of `.NormMvDataForPlot`.
#' @param x_breaks Time points at which to evaluate NAR.
#' @param y_lab Y-axis row label.
#' @return Data.frame with time, nar, event_type, event_type_label, arm.
#' @noRd
.MvOneSampleNARFrame <- function(mv_prep, x_breaks, y_lab = "Arm") {

  frames <- lapply(mv_prep$event_types, function(k) {
    proc <- .ProcessDataForPlotK(
      data = mv_prep$data,
      k = k,
      arm = NULL,
      cens_after_last = mv_prep$cens_after_last
    )
    g <- .MvNARCurveFromProc(proc)
    data.frame(
      time = x_breaks,
      nar = g(x_breaks),
      event_type = k
    )
  })
  out <- do.call(rbind, frames)
  rownames(out) <- NULL
  out <- .MvAddEventTypeLabel(
    data = out,
    event_type_labs = mv_prep$event_type_labs
  )
  out[["arm"]] <- factor(y_lab, levels = y_lab)
  out[["panel"]] <- "NAR"
  return(out)
}


#' @return ggplot or cowplot grid object.
#' @noRd
.MvPlotOneSampleMCFPanel <- function(
    mcf_df,
    color = "#C65842",
    color_lab = "Arm",
    title = NULL,
    show_auc = FALSE,
    tau = NULL,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  mcf <- time <- NULL
  q <- ggplot2::ggplot() +
    .MvGgTheme()
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
  q <- .MvApplyAxes(
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


#' Single event-type NAR panel (one sample).
#'
#' @param nar_df NAR data.frame for one event type.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis tick labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis label.
#' @return ggplot object.
#' @noRd
.MvPlotOneSampleNARPanel <- function(
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
    .MvNARTheme() +
    ggplot2::geom_text(
      ggplot2::aes(x = time, y = arm, label = nar),
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_discrete(name = NULL)
  return(q)
}


#' Single event-type MCF panel (multi-arm).
#'
#' @param mcf_df MCF data.frame for one event type.
#' @param arm_style List from `.MvResolveArmStyle`.
#' @param title Panel title.
#' @param show_auc Shade area under curves?
#' @param tau Truncation time for shading.
#' @param x_breaks X breaks.
#' @param x_lim X limits.
#' @param x_name X-axis label.
#' @param y_breaks Y breaks.
#' @param y_lim Y limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' @noRd
.MvPlotMultiArmMCFPanel <- function(
    mcf_df,
    arm_style,
    title = NULL,
    show_auc = FALSE,
    tau = NULL,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  arms <- arm_style$arms
  colors <- unname(arm_style$colors[as.character(arms)])
  labels <- unname(arm_style$color_labs[as.character(arms)])
  arm <- mcf <- time <- NULL
  q <- ggplot2::ggplot() +
    .MvGgTheme()
  if (show_auc) {
    mcf_shade <- mcf_df[mcf_df[["time"]] <= tau, , drop = FALSE]
    q <- q +
      ggplot2::geom_ribbon(
        data = mcf_shade,
        ggplot2::aes(x = time, ymin = 0, ymax = mcf, fill = arm),
        alpha = 0.5
      ) +
      ggplot2::scale_fill_manual(
        values = colors,
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
      values = colors,
      labels = labels
    )
  q <- .MvApplyAxes(
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


#' Single event-type NAR panel (multi-arm).
#'
#' @param nar_df NAR data.frame for one event type.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis tick labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis label.
#' @return ggplot object.
#' @noRd
.MvPlotMultiArmNARPanel <- function(
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
    .MvNARTheme() +
    ggplot2::geom_text(
      ggplot2::aes(x = time, y = arm, label = nar),
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_discrete(name = NULL)
  return(q)
}


#' Stack MCF and NAR panels per event type into columns.
#'
#' @param column_panels List of cowplot column grids.
#' @param title Optional overall title.
#' @return cowplot grid object.
#' @noRd
.MvAssembleEventTypeColumns <- function(column_panels, title = NULL) {

  if (length(column_panels) < 1L) {
    stop("No event-type columns to assemble.", call. = FALSE)
  }
  if (length(column_panels) == 1L) {
    out <- column_panels[[1L]]
  } else {
    out <- cowplot::plot_grid(
      plotlist = column_panels,
      nrow = 1,
      align = "h",
      axis = "tb"
    )
  }
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


#' Build combined MCF + NAR plot by event-type column.
#'
#' @param mv_prep Output of `.NormMvDataForPlot` with event type labels.
#' @param one_arm Logical; one-arm vs multi-arm layout.
#' @param arm_style Arm style list for multi-arm plots.
#' @param tau Truncation time for MCF grid and AUC shading.
#' @param x_breaks NAR evaluation times.
#' @param x_labs X-axis tick labels.
#' @param x_max X-axis upper limit.
#' @param color Line color for one-arm plots.
#' @param color_lab Legend label for one-arm plots.
#' @param show_auc Shade area under MCF curves?
#' @param scales Facet y scales.
#' @param title Plot title.
#' @param x_lim X limits.
#' @param x_name X-axis label.
#' @param y_breaks Y breaks for MCF row.
#' @param y_lim Y limits for MCF row.
#' @param y_name Y-axis label for MCF row.
#' @return ggplot or cowplot grid object.
#' @noRd
.MvBuildCombinedPlot <- function(
    mv_prep,
    one_arm,
    arm_style,
    tau,
    x_breaks,
    x_labs = NULL,
    x_max = NULL,
    color = "#C65842",
    color_lab = "Arm",
    show_auc = FALSE,
    scales = c("fixed", "free_y"),
    title = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  scales <- match.arg(scales)
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
    mcf_df_all <- MvMCFPlotFrame(
      mv_prep = mv_prep,
      eval_points = 200L,
      tau = tau
    )
    nar_df_all <- .MvOneSampleNARFrame(
      mv_prep = mv_prep,
      x_breaks = x_breaks,
      y_lab = color_lab
    )
  } else {
    mcf_df_all <- MvMultiArmMCFFrame(
      mv_prep = mv_prep,
      eval_points = 200L,
      tau = tau,
      arms = arm_style$arms
    )
    nar_df_all <- MvMultiArmNARFrame(
      mv_prep = mv_prep,
      x_breaks = x_breaks,
      arm_style = arm_style
    )
  }
  if (scales == "fixed" && is.null(y_lim)) {
    y_lim <- c(0, max(mcf_df_all[["mcf"]], na.rm = TRUE))
  }

  column_panels <- lapply(mv_prep$event_types, function(k) {
    panel_title <- .MvEventTypePanelTitle(mv_prep = mv_prep, k = k)
    mcf_df <- mcf_df_all[mcf_df_all[["event_type"]] == k, , drop = FALSE]
    nar_df <- nar_df_all[nar_df_all[["event_type"]] == k, , drop = FALSE]
    if (one_arm) {
      p_mcf <- .MvPlotOneSampleMCFPanel(
        mcf_df = mcf_df,
        color = color,
        color_lab = color_lab,
        title = panel_title,
        show_auc = show_auc,
        tau = tau,
        x_breaks = x_breaks,
        x_lim = x_lim,
        x_name = NULL,
        y_breaks = y_breaks,
        y_lim = y_lim,
        y_name = y_name
      )
      p_nar <- .MvPlotOneSampleNARPanel(
        nar_df = nar_df,
        x_breaks = x_breaks,
        x_labs = x_labs,
        x_max = x_max,
        x_name = x_name
      )
    } else {
      p_mcf <- .MvPlotMultiArmMCFPanel(
        mcf_df = mcf_df,
        arm_style = arm_style,
        title = panel_title,
        show_auc = show_auc,
        tau = tau,
        x_breaks = x_breaks,
        x_lim = x_lim,
        x_name = NULL,
        y_breaks = y_breaks,
        y_lim = y_lim,
        y_name = y_name
      )
      p_nar <- .MvPlotMultiArmNARPanel(
        nar_df = nar_df,
        x_breaks = x_breaks,
        x_labs = x_labs,
        x_max = x_max,
        x_name = x_name
      )
    }
    cowplot::plot_grid(
      p_mcf,
      p_nar,
      ncol = 1L,
      rel_heights = c(3, 1),
      align = "v",
      axis = "l"
    )
  })

  out <- .MvAssembleEventTypeColumns(
    column_panels = column_panels,
    title = title
  )
  return(out)
}


#' Facet-wrapped MCF-only multivariate plot.
#'
#' @param mv_prep Output of `.NormMvDataForPlot` with event type labels.
#' @param one_arm Logical; one-arm vs multi-arm layout.
#' @param arm_style Arm style list for multi-arm plots.
#' @param tau Truncation time.
#' @param color Line color for one-arm plots.
#' @param color_lab Legend label for one-arm plots.
#' @param show_auc Shade area under MCF curves?
#' @param scales Facet y scales.
#' @param title Plot title.
#' @param x_breaks X breaks.
#' @param x_lim X limits.
#' @param x_name X-axis label.
#' @param y_breaks Y breaks.
#' @param y_lim Y limits.
#' @param y_name Y-axis label.
#' @return ggplot object.
#' @noRd
.MvPlotMvMCFFacet <- function(
    mv_prep,
    one_arm,
    arm_style,
    tau,
    color = "#C65842",
    color_lab = "Arm",
    show_auc = FALSE,
    scales = c("fixed", "free_y"),
    title = NULL,
    x_breaks = NULL,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  scales <- match.arg(scales)
  facet_scales <- if (scales == "free_y") {
    "free_y"
  } else {
    "fixed"
  }

  if (one_arm) {
    df <- MvMCFPlotFrame(
      mv_prep = mv_prep,
      eval_points = 200L,
      tau = tau
    )
    mcf <- time <- NULL
    q <- ggplot2::ggplot() +
      .MvGgTheme()
    if (show_auc) {
      df_shade <- df[df[["time"]] <= tau, , drop = FALSE]
      q <- q +
        ggplot2::geom_ribbon(
          data = df_shade,
          ggplot2::aes(x = time, ymin = 0, ymax = mcf),
          fill = color,
          alpha = 0.5
        )
    }
    q <- q +
      ggplot2::geom_step(
        data = df,
        ggplot2::aes(x = time, y = mcf, color = "curve"),
        linewidth = 1
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = c(curve = color),
        labels = color_lab
      )
  } else {
    df <- MvMultiArmMCFFrame(
      mv_prep = mv_prep,
      eval_points = 200L,
      tau = tau,
      arms = arm_style$arms
    )
    arms <- arm_style$arms
    colors <- unname(arm_style$colors[as.character(arms)])
    labels <- unname(arm_style$color_labs[as.character(arms)])
    arm <- mcf <- time <- NULL
    q <- ggplot2::ggplot() +
      .MvGgTheme()
    if (show_auc) {
      df_shade <- df[df[["time"]] <= tau, , drop = FALSE]
      q <- q +
        ggplot2::geom_ribbon(
          data = df_shade,
          ggplot2::aes(x = time, ymin = 0, ymax = mcf, fill = arm),
          alpha = 0.5
        ) +
        ggplot2::scale_fill_manual(
          values = colors,
          guide = "none"
        )
    }
    q <- q +
      ggplot2::geom_step(
        data = df,
        ggplot2::aes(x = time, y = mcf, color = arm),
        linewidth = 1
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = colors,
        labels = labels
      )
  }

  q <- q +
    ggplot2::facet_wrap(
      ~ event_type_label,
      scales = facet_scales,
      nrow = ceiling(sqrt(mv_prep$K))
    )
  q <- .MvApplyAxes(
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


#' Apply continuous x/y scales to a ggplot.
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
.MvApplyAxes <- function(
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


#' Calculate MCF from formatted per-process data.
#'
#' @param proc Per-process data.frame.
#' @return MCF data.frame as returned by \code{CalcMCF}.
#' @noRd
.MvCalcMCFFromProc <- function(proc) {

  if (NROW(proc) < 1L || !any(proc[["status"]] == 1L)) {
    out <- data.frame(
      time = 0,
      censor = 0,
      death = 0,
      events = 0,
      nar = length(unique(proc[["idx"]])),
      haz = 0,
      surv = 1,
      event_rate = 0,
      weighted_event_rate = 0,
      mcf = 0,
      var_mcf = 0,
      se_mcf = 0
    )
    return(out)
  }

  out <- CalcMCF(
    idx = proc[["idx"]],
    status = proc[["status"]],
    time = proc[["time"]],
    jump_weights = proc[["jump_weights"]],
    calc_var = FALSE
  )
  return(out)
}


#' Build NAR stepfun from per-process formatted data.
#'
#' @param proc Per-process data.frame.
#' @return stepfun for number at risk.
#' @noRd
.MvNARCurveFromProc <- function(proc) {

  if (NROW(proc) < 1L) {
    g <- stats::stepfun(x = 0, y = c(0, 0))
    return(g)
  }

  if (!any(proc[["status"]] == 1L)) {
    g <- stats::stepfun(x = 0, y = c(0, 0))
    return(g)
  }

  fit <- CalcMCF(
    idx = proc[["idx"]],
    status = proc[["status"]],
    time = proc[["time"]],
    jump_weights = proc[["jump_weights"]],
    calc_var = FALSE
  )

  last_row <- fit[nrow(fit), , drop = FALSE]
  last_row[["time"]] <- last_row[["time"]] + 1e-8
  last_row[["nar"]] <- last_row[["nar"]] - (
    last_row[["censor"]] + last_row[["death"]]
  )
  fit <- rbind(fit, last_row)

  g <- stats::stepfun(
    x = fit[["time"]],
    y = c(length(unique(proc[["idx"]])), fit[["nar"]])
  )
  return(g)
}


# -----------------------------------------------------------------------------
# Plotting frames.
# -----------------------------------------------------------------------------

#' Multivariate MCF Plotting Frame (One Sample)
#'
#' Constructs stepped MCF values for each event type.
#'
#' @param mv_prep Output of `.NormMvDataForPlot`.
#' @param eval_points Number of evaluation points per process.
#' @param tau Truncation time for the time grid.
#' @return Data.frame with time, mcf, event_type, event_type_label.
#' @export
MvMCFPlotFrame <- function(
    mv_prep,
    eval_points = 1000,
    tau = NULL
) {

  if (is.null(tau)) {
    tau <- max(mv_prep$data[["time"]])
  }

  frames <- lapply(mv_prep$event_types, function(k) {
    proc <- .ProcessDataForPlotK(
      data = mv_prep$data,
      k = k,
      arm = NULL,
      cens_after_last = mv_prep$cens_after_last
    )
    if (NROW(proc) < 1L) {
      return(NULL)
    }
    mcf <- .MvCalcMCFFromProc(proc)
    g <- stats::stepfun(x = mcf[["time"]], y = c(0, mcf[["mcf"]]))
    times <- seq(from = 0, to = tau, length.out = eval_points)
    out <- data.frame(
      time = times,
      mcf = g(times),
      event_type = k
    )
    return(out)
  })

  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (length(frames) < 1L) {
    stop("No event processes available for plotting.", call. = FALSE)
  }

  out <- do.call(rbind, frames)
  rownames(out) <- NULL
  out <- .MvAddEventTypeLabel(
    data = out,
    event_type_labs = mv_prep$event_type_labs
  )
  return(out)
}


#' Multivariate multi-arm MCF plotting frame.
#'
#' @param mv_prep Output of `.NormMvDataForPlot`.
#' @param eval_points Number of evaluation points per arm and process.
#' @param tau Truncation time for the time grid.
#' @param arms Sorted vector of arm codes to include.
#' @return Data.frame with time, mcf, arm, event_type, event_type_label.
#' @noRd
MvMultiArmMCFFrame <- function(
    mv_prep,
    eval_points = 200,
    tau = NULL,
    arms = NULL
) {

  if (is.null(tau)) {
    tau <- max(mv_prep$data[["time"]])
  }
  if (is.null(arms)) {
    arms <- sort(unique(mv_prep$data[["arm"]]))
  }

  times <- seq(from = 0, to = tau, length.out = eval_points)
  frames <- list()

  for (k in mv_prep$event_types) {
    for (a in arms) {
      proc <- .ProcessDataForPlotK(
        data = mv_prep$data,
        k = k,
        arm = a,
        cens_after_last = mv_prep$cens_after_last
      )
      if (NROW(proc) < 1L) {
        next
      }
      mcf <- .MvCalcMCFFromProc(proc)
      g <- stats::stepfun(x = mcf[["time"]], y = c(0, mcf[["mcf"]]))
      frames[[length(frames) + 1L]] <- data.frame(
        time = times,
        mcf = g(times),
        arm = a,
        event_type = k
      )
    }
  }

  if (length(frames) < 1L) {
    stop("No event processes available for plotting.", call. = FALSE)
  }

  out <- do.call(rbind, frames)
  rownames(out) <- NULL
  out[["arm"]] <- factor(out[["arm"]], levels = arms)
  out <- .MvAddEventTypeLabel(
    data = out,
    event_type_labs = mv_prep$event_type_labs
  )
  return(out)
}


#' Multivariate multi-arm NAR plotting frame.
#'
#' @param mv_prep Output of `.NormMvDataForPlot`.
#' @param x_breaks Time points at which to evaluate NAR.
#' @param arm_style List from `.MvResolveArmStyle`.
#' @return Data.frame with time, nar, arm, event_type, event_type_label.
#' @noRd
MvMultiArmNARFrame <- function(
    mv_prep,
    x_breaks,
    arm_style
) {

  arms <- arm_style$arms
  lab_vec <- unname(arm_style$color_labs[as.character(arms)])
  frames <- list()

  for (k in mv_prep$event_types) {
    for (a in arms) {
      proc <- .ProcessDataForPlotK(
        data = mv_prep$data,
        k = k,
        arm = a,
        cens_after_last = mv_prep$cens_after_last
      )
      g <- .MvNARCurveFromProc(proc)
      frames[[length(frames) + 1L]] <- data.frame(
        time = x_breaks,
        nar = g(x_breaks),
        arm = a,
        event_type = k
      )
    }
  }

  if (length(frames) < 1L) {
    stop("No event processes available for NAR plotting.", call. = FALSE)
  }

  out <- do.call(rbind, frames)
  rownames(out) <- NULL
  out[["arm"]] <- factor(
    out[["arm"]],
    levels = arms,
    labels = lab_vec
  )
  out <- .MvAddEventTypeLabel(
    data = out,
    event_type_labs = mv_prep$event_type_labs
  )
  return(out)
}


# -----------------------------------------------------------------------------
# Curve constructors.
# -----------------------------------------------------------------------------

#' Multivariate MCF Curves
#'
#' Return stepfun objects evaluating the mean cumulative function at each
#' event type for one sample.
#'
#' @param data Data.frame in long format with an \code{event_type} column.
#' @param cens_after_last Censor subjects without a terminal row after their
#'   last recurrent event?
#' @param event_type_name Name of event type column.
#' @param idx_name Name of subject index column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Named list of \code{stepfun} objects, one per event type.
#' @export
MvMCFCurve <- function(
    data,
    cens_after_last = TRUE,
    event_type_name = "event_type",
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
) {

  mv_prep <- .NormMvDataForPlot(
    data = data,
    arm_name = NULL,
    cens_after_last = cens_after_last,
    event_type_name = event_type_name,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name
  )

  curves <- lapply(mv_prep$event_types, function(k) {
    proc <- .ProcessDataForPlotK(
      data = mv_prep$data,
      k = k,
      arm = NULL,
      cens_after_last = mv_prep$cens_after_last
    )
    mcf <- .MvCalcMCFFromProc(proc)
    g <- stats::stepfun(x = mcf[["time"]], y = c(0, mcf[["mcf"]]))
    return(g)
  })
  names(curves) <- as.character(mv_prep$event_types)
  return(curves)
}


#' Multivariate Number at Risk Curves
#'
#' @param data Data.frame in long format with an \code{event_type} column.
#' @param cens_after_last Censor subjects without a terminal row after their
#'   last recurrent event?
#' @param event_type_name Name of event type column.
#' @param idx_name Name of subject index column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Named list of \code{stepfun} objects, one per event type.
#' @export
MvNARCurve <- function(
    data,
    cens_after_last = TRUE,
    event_type_name = "event_type",
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
) {

  mv_prep <- .NormMvDataForPlot(
    data = data,
    arm_name = NULL,
    cens_after_last = cens_after_last,
    event_type_name = event_type_name,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name
  )

  curves <- lapply(mv_prep$event_types, function(k) {
    proc <- .ProcessDataForPlotK(
      data = mv_prep$data,
      k = k,
      arm = NULL,
      cens_after_last = mv_prep$cens_after_last
    )
    g <- .MvNARCurveFromProc(proc)
    return(g)
  })
  names(curves) <- as.character(mv_prep$event_types)
  return(curves)
}



# -----------------------------------------------------------------------------
# Multivariate plots.
# -----------------------------------------------------------------------------

#' Plot Multivariate Mean Cumulative Functions
#'
#' Plot mean cumulative functions by treatment arm and event type. Supports
#' one to \eqn{K} arms, optional numbers at risk below each event-type column
#' (\code{show_nar = TRUE}), and optional area shading (\code{show_auc = TRUE}).
#' Curves are empirical per-process MCFs (the same objects underlying
#' \code{@MCF} from \code{\link{CompareMvAUCs}}).
#'
#' @param data Data.frame in long format with \code{event_type}.
#' @param which_arm Optional arm code to plot a single arm; \code{NULL} plots
#'   all arms present.
#' @param arm_name Name of arm column.
#' @param color_labs Legend labels for multi-arm plots (length = number of
#'   arms; default \code{Arm <code>}).
#' @param colors Line and fill colors for multi-arm plots (length >= number of
#'   arms).
#' @param color Line color when a single arm is plotted.
#' @param color_lab Legend label when a single arm is plotted.
#' @param cens_after_last Censor subjects without a terminal row after their
#'   last recurrent event?
#' @param event_type_name Name of event type column.
#' @param event_type_labs Optional labels for each event type.
#' @param idx_name Name of subject index column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param tau Truncation time for the time grid.
#' @param title Plot title.
#' @param scales Facet y scales when \code{show_nar = FALSE}:
#'   \code{"fixed"} or \code{"free_y"}.
#' @param x_breaks X-axis breaks for the NAR row (derived from \code{tau} when
#'   \code{NULL} and \code{show_nar = TRUE}).
#' @param x_labs X-axis tick labels for the NAR row.
#' @param x_max X-axis upper limit for the NAR row.
#' @param show_nar Include numbers at risk below each event-type MCF panel?
#' @param show_auc Shade the area under each MCF curve from 0 to \code{tau}?
#' @param x_lim X-axis limits.
#' @param x_name X-axis label.
#' @param y_breaks Y-axis breaks.
#' @param y_lim Y-axis limits.
#' @param y_name Y-axis label.
#' @return A \code{ggplot} object, or a cowplot layout when \code{show_nar =
#'   TRUE}.
#' @export
PlotMvMCFs <- function(
    data,
    which_arm = NULL,
    arm_name = "arm",
    color_labs = NULL,
    colors = NULL,
    color = "#C65842",
    color_lab = "Arm",
    cens_after_last = TRUE,
    event_type_name = "event_type",
    event_type_labs = NULL,
    idx_name = "idx",
    status_name = "status",
    time_name = "time",
    tau = NULL,
    title = NULL,
    scales = c("fixed", "free_y"),
    x_breaks = NULL,
    x_labs = NULL,
    x_max = NULL,
    show_nar = TRUE,
    show_auc = FALSE,
    x_lim = NULL,
    x_name = "Time",
    y_breaks = NULL,
    y_lim = NULL,
    y_name = "Mean Cumulative Count"
) {

  scales <- match.arg(scales)
  data <- .MvSubsetWhichArm(
    data = data,
    which_arm = which_arm,
    arm_name = arm_name
  )
  arms <- .MvUniqueArms(data = data, arm_name = arm_name)
  n_arms <- length(arms)
  one_arm <- n_arms <= 1L
  arm_style <- NULL
  if (!one_arm) {
    arm_style <- .MvResolveArmStyle(
      arms = arms,
      color_labs = color_labs,
      colors = colors
    )
  }

  if (!event_type_name %in% names(data)) {
    if (n_arms > 2L) {
      stop(
        "Multivariate data (event_type column) is required for K > 2 arm plots.",
        call. = FALSE
      )
    }
    if (one_arm) {
      out <- PlotMCFs(
        data = data,
        color = color,
        color_lab = color_lab,
        idx_name = idx_name,
        status_name = status_name,
        time_name = time_name,
        title = title,
        tau = tau,
        x_breaks = x_breaks,
        x_lim = x_lim,
        x_name = x_name,
        y_breaks = y_breaks,
        y_lim = y_lim,
        y_name = y_name,
        show_nar = show_nar,
        show_auc = show_auc
      )
      return(out)
    }
    if (show_auc) {
      out <- PlotMCFs(
        data = data,
        which_arm = arms[1L],
        arm_name = arm_name,
        color = arm_style$colors[[as.character(arms[1L])]],
        color_lab = arm_style$color_labs[[as.character(arms[1L])]],
        idx_name = idx_name,
        status_name = status_name,
        time_name = time_name,
        title = title,
        tau = tau,
        x_breaks = x_breaks,
        x_lim = x_lim,
        x_name = x_name,
        y_breaks = y_breaks,
        y_lim = y_lim,
        y_name = y_name,
        show_nar = show_nar,
        show_auc = TRUE
      )
      return(out)
    }
    out <- PlotMCFs(
      data = data,
      arm_name = arm_name,
      color_labs = unname(arm_style$color_labs),
      colors = unname(arm_style$colors),
      idx_name = idx_name,
      status_name = status_name,
      tau = tau,
      time_name = time_name,
      title = title,
      x_breaks = x_breaks,
      x_lim = x_lim,
      x_name = x_name,
      y_breaks = y_breaks,
      y_lim = y_lim,
      y_name = y_name,
      show_nar = show_nar,
      show_auc = FALSE
    )
    return(out)
  }

  mv_prep <- .NormMvDataForPlot(
    data = data,
    arm_name = arm_name,
    cens_after_last = cens_after_last,
    event_type_name = event_type_name,
    idx_name = idx_name,
    status_name = status_name,
    time_name = time_name
  )
  mv_prep <- .MvPrepEventTypeLabs(
    mv_prep = mv_prep,
    event_type_labs = event_type_labs
  )

  arms <- .MvUniqueArms(data = mv_prep$data, arm_name = "arm")
  n_arms <- length(arms)
  one_arm <- n_arms <= 1L
  if (!one_arm) {
    arm_style <- .MvResolveArmStyle(
      arms = arms,
      color_labs = color_labs,
      colors = colors
    )
  }

  if (is.null(x_lim) || is.null(x_lim[2])) {
    x_max_data <- max(mv_prep$data[["time"]])
  } else {
    x_max_data <- x_lim[2]
  }
  if (is.null(x_max)) {
    x_max <- x_max_data
  }
  if (is.null(tau)) {
    tau <- x_max_data
  }

  if (show_nar) {
    if (is.null(x_breaks)) {
      x_breaks <- .MvDefaultXBreaks(tau)
    }
    out <- .MvBuildCombinedPlot(
      mv_prep = mv_prep,
      one_arm = one_arm,
      arm_style = arm_style,
      tau = tau,
      x_breaks = x_breaks,
      x_labs = x_labs,
      x_max = x_max,
      color = color,
      color_lab = color_lab,
      show_auc = show_auc,
      scales = scales,
      title = title,
      x_lim = x_lim,
      x_name = x_name,
      y_breaks = y_breaks,
      y_lim = y_lim,
      y_name = y_name
    )
    return(out)
  }

  out <- .MvPlotMvMCFFacet(
    mv_prep = mv_prep,
    one_arm = one_arm,
    arm_style = arm_style,
    tau = tau,
    color = color,
    color_lab = color_lab,
    show_auc = show_auc,
    scales = scales,
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
