# Purpose: Format multivariate recurrent events data.
# Updated: 2026-06-16

#' Format Multivariate Data
#'
#' Prepare long-format multivariate recurrent events data for analysis.
#'
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param cens_after_last Censor after last event if no terminal row?
#' @param covars Optional subject-level covariates.
#' @param event_type_name Name of event type column.
#' @param idx_name Name of subject index column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return List with formatted `data`, integer vector `event_types`, and
#'   integer `K`.
#' @noRd
FormatMvData <- function(
    data,
    arm_name = "arm",
    cens_after_last = TRUE,
    covars = NULL,
    event_type_name = "event_type",
    idx_name = "idx",
    status_name = "status",
    time_name = "time"
) {

  # Rename core columns.
  key_cols <- c(arm_name, idx_name, status_name, time_name)
  has_event_type <- event_type_name %in% names(data)
  if (has_event_type) {
    key_cols <- c(key_cols, event_type_name)
  }

  data <- data %>%
    dplyr::select(dplyr::all_of(key_cols)) %>%
    dplyr::rename(
      arm = {{ arm_name }},
      idx = {{ idx_name }},
      status = {{ status_name }},
      time = {{ time_name }}
    )

  if (has_event_type) {
    data <- data %>% dplyr::rename(event_type = {{ event_type_name }})
  } else {
    data$event_type <- 1L
    has_event_type <- FALSE
  }

  # Univariate shortcut: delegate to FormatData when no explicit event types.
  if (!has_event_type) {
    data <- FormatData(
      data = data,
      arm_name = "arm",
      cens_after_last = cens_after_last,
      covars = covars,
      idx_name = "idx",
      status_name = "status",
      time_name = "time"
    )
    data$event_type <- 1L
    out <- list(
      data = data,
      event_types = 1L,
      K = 1L,
      univariate_path = TRUE
    )
    return(out)
  }

  # Multivariate path.
  data <- ConvertIdxToInt(data)

  if (!is.null(covars)) {
    data <- cbind(data, covars)
  }

  data$jump_weights <- 1

  # Recurrent rows must have event_type; terminals may be NA (global).
  recurrent <- data[["status"]] == 1
  if (any(is.na(data[["event_type"]][recurrent]))) {
    stop("Recurrent events (status == 1) must have a non-missing event_type.")
  }

  typed_k <- data[["event_type"]][!is.na(data[["event_type"]])]
  event_types <- sort(unique(as.integer(typed_k)))
  K <- length(event_types)

  out <- list(
    data = data,
    event_types = event_types,
    K = K,
    univariate_path = FALSE,
    cens_after_last = cens_after_last
  )

  .CheckMvTerminalEvents(
    data = data,
    event_types = event_types,
    cens_after_last = cens_after_last
  )

  return(out)
}


#' Subject indices at risk for one event process.
#'
#' A subject is at risk for process \eqn{k} if they have at least one row with
#' non-missing \code{event_type = k} (recurrent or typed terminal).
#'
#' @param data Formatted multivariate data.frame.
#' @param k Event type.
#' @return Integer vector of subject indices.
#' @noRd
.MvAtRiskIdx <- function(data, k) {
  idx_col <- data[["idx"]]
  et_col <- data[["event_type"]]
  at_risk <- !is.na(et_col) & et_col == k
  out <- unique(idx_col[at_risk])
  return(out)
}


#' Identify the observation-terminating row for one event process.
#'
#' @param subj Subject-level data.frame.
#' @param k Event type.
#' @return List with `terminal` (0 or 1 row) and `recurrent` data.frames.
#' @noRd
.TerminalForProcessK <- function(subj, k) {

  is_recurrent <- subj[["status"]] == 1L
  is_recurrent[is.na(is_recurrent)] <- FALSE
  match_type <- subj[["event_type"]] == k
  match_type[is.na(match_type)] <- FALSE
  recurrent <- subj[is_recurrent & match_type, , drop = FALSE]

  is_terminal <- subj[["status"]] %in% c(0L, 2L)
  typed_idx <- is_terminal & (subj[["event_type"]] == k)
  typed_idx[is.na(typed_idx)] <- FALSE
  typed_term <- subj[typed_idx, , drop = FALSE]

  global_idx <- is_terminal & is.na(subj[["event_type"]])
  global_term <- subj[global_idx, , drop = FALSE]

  if (NROW(typed_term) > 1) {
    sid <- subj[["idx"]][1]
    stop(
      "Subject ", sid, " has multiple observation terminating events ",
      "for event type ", k, ".",
      call. = FALSE
    )
  }
  if (NROW(global_term) > 1) {
    sid <- subj[["idx"]][1]
    stop(
      "Subject ", sid, " has multiple global observation terminating events.",
      call. = FALSE
    )
  }

  if (NROW(typed_term) > 0) {
    terminal <- typed_term
  } else if (NROW(global_term) > 0) {
    terminal <- global_term
  } else {
    terminal <- NULL
  }

  if (!is.null(terminal)) {
    term_time <- terminal[["time"]][1]
    keep_rec <- recurrent[["time"]] <= term_time
    keep_rec[is.na(keep_rec)] <- FALSE
    recurrent <- recurrent[keep_rec, , drop = FALSE]
  }

  return(list(
    terminal = terminal,
    recurrent = recurrent
  ))
}


#' Validate observation-terminating events per at-risk event type.
#'
#' @param data Formatted multivariate data.
#' @param event_types Event type levels.
#' @param cens_after_last Add censoring after the last recurrent event?
#' @return Invisible NULL.
#' @noRd
.CheckMvTerminalEvents <- function(
    data,
    event_types,
    cens_after_last
) {

  missing_term <- character(0)
  missing_any <- character(0)

  for (k in event_types) {
    at_risk_idx <- .MvAtRiskIdx(data = data, k = k)
    for (sid in at_risk_idx) {
      subj <- data[data[["idx"]] == sid, , drop = FALSE]
      parts <- .TerminalForProcessK(subj = subj, k = k)
      has_rec <- NROW(parts$recurrent) > 0
      has_term <- !is.null(parts$terminal)
      label <- paste0("idx ", sid, ", event_type ", k)

      if (has_term) {
        next
      }
      if (has_rec && !cens_after_last) {
        missing_term <- c(missing_term, label)
      } else if (!has_rec) {
        missing_any <- c(missing_any, label)
      }
    }
  }

  if (length(missing_term) > 0 && !cens_after_last) {
    warning(
      "Patients without observation terminating events were found for: ",
      paste(missing_term, collapse = "; "),
      ".",
      call. = FALSE
    )
  }

  if (length(missing_any) > 0) {
    warning(
      "At-risk patients with no recurrent events and no observation ",
      "terminating event were found for: ",
      paste(missing_any, collapse = "; "),
      ".",
      call. = FALSE
    )
  }

  return(invisible(NULL))
}


#' Subset and Format Data for a Single Event Process
#'
#' @param data Formatted multivariate data.frame.
#' @param k Event type.
#' @param cens_after_last Censor after last event if no terminal row?
#' @return Per-process data.frame with columns idx, status, time, jump_weights.
#' @noRd
.SubsetProcessK <- function(
    data,
    k,
    cens_after_last = TRUE
) {

  empty_proc <- function() {
    return(data.frame(
      idx = integer(),
      status = integer(),
      time = numeric(),
      jump_weights = numeric()
    ))
  }

  if (NROW(data) < 1) {
    return(empty_proc())
  }

  subjects <- unique(data[["idx"]])
  if (length(subjects) < 1) {
    return(empty_proc())
  }

  proc_rows <- list()

  for (sid in subjects) {
    subj <- data[data[["idx"]] == sid, , drop = FALSE]
    if (NROW(subj) < 1) {
      next
    }

    parts <- .TerminalForProcessK(subj = subj, k = k)
    recurrent <- parts$recurrent
    terminal <- parts$terminal

    if (!is.null(terminal)) {
      subj_proc <- rbind(recurrent, terminal)
    } else {
      subj_proc <- recurrent
    }

    if (NROW(subj_proc) > 0) {
      proc_rows[[as.character(sid)]] <- subj_proc
    }
  }

  if (length(proc_rows) < 1) {
    return(empty_proc())
  }

  proc <- do.call(rbind, proc_rows)
  rownames(proc) <- NULL

  keep_cols <- intersect(
    c("idx", "status", "time", "jump_weights", "orig_idx"),
    names(proc)
  )
  proc <- proc[, keep_cols, drop = FALSE]
  proc[["obs_end"]] <- 1L * (proc[["status"]] != 1L)
  ord <- order(proc[["idx"]], proc[["time"]], proc[["obs_end"]])
  proc <- proc[ord, , drop = FALSE]

  split_data <- split(x = proc, f = proc[["idx"]])
  format_data <- lapply(split_data, function(x) {
    FormatSubj(x, cens_after_last = cens_after_last)
  })
  proc_out <- do.call(rbind, format_data)

  if (is.null(proc_out) || NROW(proc_out) < 1) {
    return(empty_proc())
  }

  rownames(proc_out) <- NULL
  if ("obs_end" %in% names(proc_out)) {
    proc_out[["obs_end"]] <- NULL
  }

  keep_idx <- proc_out %>%
    dplyr::group_by(.data[["idx"]]) %>%
    dplyr::summarise(
      n_term = sum(.data[["status"]] %in% c(0L, 2L)),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data[["n_term"]] == 1L) %>%
    dplyr::pull(.data[["idx"]])
  proc_out <- proc_out[proc_out[["idx"]] %in% keep_idx, , drop = FALSE]
  rownames(proc_out) <- NULL

  if (NROW(proc_out) < 1) {
    return(empty_proc())
  }

  return(proc_out)
}
