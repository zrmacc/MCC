# Purpose: Bootstrap inference for multivariate AUC estimator.
# Updated: 2026-06-16

#' Bootstrap Confidence Intervals for Multivariate Contrasts
#'
#' @param sim Bootstrap simulation results.
#' @param obs_stats Observed contrasts.
#' @param alpha Type I error.
#' @return Data.frame of bootstrap confidence intervals.
#' @noRd
BootCIsMv <- function(
    sim,
    obs_stats,
    alpha = 0.05
) {

  alpha2 <- alpha / 2
  event_types <- unique(obs_stats$event_type)
  out <- NULL

  for (et in event_types) {
    col_nm <- paste0("boot_diff_", et)
    if (!col_nm %in% names(sim)) {
      next
    }
    ci_diff <- stats::quantile(
      x = sim[[col_nm]],
      probs = c(alpha2, 1 - alpha2)
    )
    names(ci_diff) <- NULL
    obs_row <- obs_stats[obs_stats$event_type == et, , drop = FALSE]
    row <- data.frame(
      method = "bootstrap",
      contrast = "A1-A0",
      event_type = et,
      observed = obs_row$observed[1],
      se = stats::sd(sim[[col_nm]]),
      lower = ci_diff[1],
      upper = ci_diff[2],
      stringsAsFactors = FALSE
    )
    out <- rbind(out, row)
  }

  return(out)
}


#' Assemble Multivariate CIs and P-values with Optional Bootstrap
#'
#' @param obs_stats Observed contrast data.frame.
#' @param data Formatted data.
#' @param event_types Event types.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Bootstrap replicates.
#' @param cens_after_last Censor after last event?
#' @param covar_cols Covariate columns or NULL.
#' @param process_weights Optional contrast weights.
#' @return List with cis, pvals, sim_reps.
#' @noRd
.AddMvResamplingResults <- function(
    obs_stats,
    data,
    event_types,
    tau,
    alpha,
    reps,
    cens_after_last,
    covar_cols,
    process_weights
) {

  contrast <- observed <- p <- event_type <- NULL
  cis <- obs_stats %>% dplyr::select(-dplyr::any_of("p"))
  cis <- data.frame(method = "asymptotic", cis)
  pvals <- obs_stats %>% dplyr::select(dplyr::any_of(c(
    "contrast", "event_type", "observed", "p"
  )))
  pvals <- data.frame(method = "asymptotic", pvals)
  sim_reps <- list()

  if (!is.null(reps)) {
    boot_sim <- BootSimMvAug(
      data = data,
      event_types = event_types,
      obs_stats = obs_stats,
      tau = tau,
      alpha = alpha,
      reps = reps,
      cens_after_last = cens_after_last,
      covar_cols = covar_cols,
      process_weights = process_weights
    )
    sim_reps$boot_sim <- boot_sim
    boot_cis <- BootCIsMv(
      sim = boot_sim,
      obs_stats = obs_stats,
      alpha = alpha
    )
    cis <- rbind(cis, boot_cis)
    cis <- cis[order(cis$event_type, cis$method), ]

    event_types <- unique(obs_stats$event_type)
    boot_pvals <- NULL
    for (et in event_types) {
      sign_col <- paste0("is_diff_sign_", et)
      obs_row <- obs_stats[obs_stats$event_type == et, , drop = FALSE]
      boot_pvals <- rbind(
        boot_pvals,
        data.frame(
          method = "bootstrap",
          contrast = "A1-A0",
          event_type = et,
          observed = obs_row$observed[1],
          p = CalcP(boot_sim[[sign_col]]),
          stringsAsFactors = FALSE
        )
      )
    }
    pvals <- rbind(pvals, boot_pvals)
  }

  out <- list(
    cis = cis,
    pvals = pvals,
    sim_reps = sim_reps
  )
  return(out)
}


#' Bootstrap Inference for Multivariate Augmentation Estimator
#'
#' @param data Formatted data.
#' @param event_types Event types.
#' @param obs_stats Observed contrasts.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param reps Bootstrap replicates.
#' @param cens_after_last Censor after last event?
#' @param covar_cols Covariate columns or NULL.
#' @param process_weights Optional contrast weights.
#' @return Data.frame of bootstrap replicates.
#' @noRd
BootSimMvAug <- function(
    data,
    event_types,
    obs_stats,
    tau,
    alpha = 0.05,
    reps = 2000,
    cens_after_last = TRUE,
    covar_cols = NULL,
    process_weights = NULL
) {

  data0 <- data[data[["arm"]] == 0, , drop = FALSE]
  data1 <- data[data[["arm"]] == 1, , drop = FALSE]
  K <- length(event_types)

  Loop <- function(b) {
    boot0 <- GroupBoot(data0)
    boot1 <- GroupBoot(data1)
    boot <- rbind(boot1, boot0)
    boot_stats <- CalcMvAugAUC(
      data = boot,
      event_types = event_types,
      tau = tau,
      alpha = alpha,
      cens_after_last = cens_after_last,
      covar_cols = covar_cols,
      process_weights = process_weights,
      return_areas = FALSE
    )
    out <- as.numeric(boot_stats$contrasts$observed)
    names(out) <- paste0("boot_diff_", event_types)
    for (k in seq_len(K)) {
      et <- event_types[k]
      obs_k <- obs_stats$observed[obs_stats$event_type == et][1]
      boot_k <- boot_stats$contrasts$observed[boot_stats$contrasts$event_type == et][1]
      out[paste0("is_diff_sign_", et)] <- sign(boot_k) != sign(obs_k)
    }
    return(out)
  }

  sim <- lapply(seq_len(reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  return(sim)
}
