# Test helpers for multivariate eligibility simulation DGP utilities.
# Mirrors Multivariate/scripts/utils/multivariate_eligibility_dgp.R for use
# when the Multivariate subproject is not on the search path.

#' All eligibility scenario names.
#'
#' @return Character vector of length 4.
#' @noRd
MvEligScenarioNames <- function() {
  out <- c(
    "elig_no_covar",
    "elig_covar_event",
    "elig_covar_death",
    "elig_covar_types"
  )
  return(out)
}


#' Split integer n into three parts summing to n (deterministic).
#'
#' @param n Integer count to split.
#' @return Integer vector of length 3.
#' @noRd
SplitThirds <- function(n) {
  if (n < 3L) {
    stop("n must be at least 3 to assign eligibility thirds.", call. = FALSE)
  }
  base <- n %/% 3L
  rem <- n %% 3L
  out <- rep(base, 3L)
  if (rem >= 1L) {
    out[1L] <- out[1L] + 1L
  }
  if (rem >= 2L) {
    out[2L] <- out[2L] + 1L
  }
  return(out)
}


#' Assign eligibility patterns within each arm.
#'
#' @param n_per_arm Subjects per arm.
#' @param fractions Target fractions for both, k1_only, k2_only (sum to 1).
#' @return Data.frame with columns idx, arm, elig_pattern.
#' @noRd
AssignMvEligibility <- function(
    n_per_arm,
    fractions = c(1 / 3, 1 / 3, 1 / 3)
) {
  if (length(fractions) != 3L) {
    stop("fractions must have length 3.", call. = FALSE)
  }
  patterns <- c("both", "k1_only", "k2_only")
  arms <- c(0L, 1L)
  out_list <- vector("list", length(arms))
  next_idx <- 1L
  for (a in seq_along(arms)) {
    counts <- SplitThirds(n_per_arm)
    arm_idx <- next_idx:(next_idx + n_per_arm - 1L)
    next_idx <- next_idx + n_per_arm
    pat_vec <- rep(patterns, times = counts)
    out_list[[a]] <- data.frame(
      idx = arm_idx,
      arm = arms[a],
      elig_pattern = pat_vec,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  return(out)
}


#' Eligible event types for an eligibility pattern name.
#'
#' @param pattern One of both, k1_only, k2_only.
#' @return Integer vector of event type indices.
#' @noRd
.EligibleTypesForPattern <- function(pattern) {
  if (pattern == "both") {
    return(c(1L, 2L))
  }
  if (pattern == "k1_only") {
    return(1L)
  }
  if (pattern == "k2_only") {
    return(2L)
  }
  stop("Unknown eligibility pattern: ", pattern, call. = FALSE)
}


#' Post-process simulated data to enforce event-type eligibility.
#'
#' @param data Long-format data from \code{GenData}.
#' @param elig_map Data.frame with idx and elig_pattern columns.
#' @return Long-format data with ineligible recurrent rows removed.
#' @noRd
ApplyMvEventEligibility <- function(data, elig_map) {
  if (!all(c("idx", "elig_pattern") %in% names(elig_map))) {
    stop("elig_map must contain idx and elig_pattern.", call. = FALSE)
  }
  idx_col <- data[["idx"]]
  subj_ids <- unique(idx_col)
  out_rows <- vector("list", length(subj_ids))
  for (i in seq_along(subj_ids)) {
    sid <- subj_ids[i]
    subj <- data[idx_col == sid, , drop = FALSE]
    pat <- elig_map$elig_pattern[elig_map$idx == sid][1L]
    if (is.na(pat)) {
      stop("Missing eligibility pattern for idx ", sid, ".", call. = FALSE)
    }
    eligible <- .EligibleTypesForPattern(pat)
    is_global <- subj[["status"]] %in% c(0L, 2L) & is.na(subj[["event_type"]])
    global_rows <- subj[is_global, , drop = FALSE]
    if (NROW(global_rows) != 1L) {
      stop("Subject ", sid, " must have exactly one global terminal row.", call. = FALSE)
    }
    global_time <- global_rows[["time"]][1L]
    global_status <- global_rows[["status"]][1L]
    is_rec <- subj[["status"]] == 1L
    is_rec[is.na(is_rec)] <- FALSE
    rec <- subj[is_rec, , drop = FALSE]
    if (NROW(rec) > 0L) {
      keep_rec <- rec[["event_type"]] %in% eligible
      rec <- rec[keep_rec, , drop = FALSE]
    }
    meta_cols <- setdiff(names(subj), c("time", "status", "event_type"))
    typed_terms <- vector("list", length(eligible))
    for (k in seq_along(eligible)) {
      et <- eligible[k]
      has_k <- any(subj[["event_type"]] == et, na.rm = TRUE)
      if (!has_k) {
        typed_row <- global_rows[1L, meta_cols, drop = FALSE]
        typed_row$time <- global_time
        typed_row$status <- global_status
        typed_row$event_type <- et
        typed_terms[[k]] <- typed_row
      }
    }
    typed_terms <- typed_terms[!vapply(typed_terms, is.null, logical(1L))]
    subj_out <- rbind(rec, do.call(rbind, typed_terms), global_rows)
    rownames(subj_out) <- NULL
    subj_out$elig_pattern <- pat
    out_rows[[i]] <- subj_out
  }
  out <- do.call(rbind, out_rows)
  rownames(out) <- NULL
  out <- out[order(out$idx, out$time, out$event_type, na.last = TRUE), ]
  rownames(out) <- NULL
  return(out)
}


#' Build covariate design for an eligibility scenario.
#'
#' @param n_per_arm Subjects per arm.
#' @param scenario_name Eligibility scenario name.
#' @return Data.frame with intercept, arm, and optional covariates.
#' @noRd
MakeMvEligDesign <- function(n_per_arm, scenario_name) {
  if (!scenario_name %in% MvEligScenarioNames()) {
    stop("Unknown eligibility scenario: ", scenario_name, call. = FALSE)
  }
  elig <- AssignMvEligibility(n_per_arm = n_per_arm)
  n <- nrow(elig)
  out <- data.frame(
    idx = elig$idx,
    intercept = 1,
    arm = elig$arm,
    elig_pattern = elig$elig_pattern,
    stringsAsFactors = FALSE
  )
  if (scenario_name == "elig_covar_event" || scenario_name == "elig_covar_death") {
    out$X <- stats::rnorm(n, mean = 0, sd = 1)
  } else if (scenario_name == "elig_covar_types") {
    out$X1 <- stats::rnorm(n, mean = 0, sd = 1)
    out$X2 <- stats::rnorm(n, mean = 0, sd = 1)
  }
  return(out)
}


#' Scenario configuration for eligibility simulations.
#'
#' @param scenario_name Eligibility scenario name.
#' @param study Either t1e or bias.
#' @return List with GenData arguments and augmentation metadata.
#' @noRd
.MvEligScenarioConfig <- function(scenario_name, study = c("t1e", "bias")) {
  study <- match.arg(study)
  if (!scenario_name %in% MvEligScenarioNames()) {
    stop("Unknown eligibility scenario: ", scenario_name, call. = FALSE)
  }
  log_075 <- log(0.75)
  log_125 <- log(1.25)
  arm_event <- if (study == "bias") {
    log_075
  } else {
    0
  }
  shell <- list(
    base_event_rate = c(1.0, 1.0),
    event_frailty_variance = 1.0,
    death_frailty_variance = 1.0,
    scenario = scenario_name,
    study = study
  )
  if (scenario_name == "elig_no_covar") {
    shell$beta_event <- c(0, arm_event)
    shell$beta_death <- c(0, 0)
    shell$use_augmentation <- FALSE
    shell$covar_cols <- character(0)
  } else if (scenario_name == "elig_covar_event") {
    shell$beta_event <- c(0, arm_event, log_125)
    shell$beta_death <- c(0, 0, 0)
    shell$use_augmentation <- TRUE
    shell$covar_cols <- "X"
  } else if (scenario_name == "elig_covar_death") {
    shell$beta_event <- c(0, arm_event, 0)
    shell$beta_death <- c(0, 0, log_125)
    shell$use_augmentation <- TRUE
    shell$covar_cols <- "X"
  } else if (scenario_name == "elig_covar_types") {
    beta_event <- matrix(0, nrow = 4L, ncol = 2L)
    colnames(beta_event) <- c("k1", "k2")
    rownames(beta_event) <- c("intercept", "arm", "X1", "X2")
    beta_event["arm", ] <- arm_event
    beta_event["X1", "k1"] <- log_125
    beta_event["X2", "k2"] <- log_075
    shell$beta_event <- beta_event
    shell$beta_death <- c(0, 0, 0, 0)
    shell$use_augmentation <- TRUE
    shell$covar_cols <- c("X1", "X2")
  }
  return(shell)
}


#' Simulate multivariate data under an eligibility DGP.
#'
#' @param scenario_name Eligibility scenario name.
#' @param study Either t1e or bias.
#' @param n_per_arm Subjects per arm.
#' @param tau_gen Administrative truncation for data generation.
#' @return Long-format data.frame with arm and elig_pattern columns.
#' @noRd
SimMvEligData <- function(
    scenario_name,
    study = c("t1e", "bias"),
    n_per_arm,
    tau_gen = 4
) {
  study <- match.arg(study)
  cfg <- .MvEligScenarioConfig(
    scenario_name = scenario_name,
    study = study
  )
  design <- MakeMvEligDesign(
    n_per_arm = n_per_arm,
    scenario_name = scenario_name
  )
  elig_map <- design[, c("idx", "arm", "elig_pattern"), drop = FALSE]
  cov_cols <- setdiff(names(design), c("idx", "elig_pattern"))
  data <- GenData(
    covariates = design[, cov_cols, drop = FALSE],
    n_event_types = 2L,
    base_event_rate = cfg$base_event_rate,
    beta_event = cfg$beta_event,
    beta_death = cfg$beta_death,
    event_frailty_variance = cfg$event_frailty_variance,
    death_frailty_variance = cfg$death_frailty_variance,
    tau = tau_gen
  )
  if (!"arm" %in% names(data)) {
    stop("Simulated data are missing arm column.", call. = FALSE)
  }
  data <- ApplyMvEventEligibility(
    data = data,
    elig_map = elig_map
  )
  data$eligibility_design <- "thirds"
  return(data)
}
