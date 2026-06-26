# Purpose: Multivariate AUMCF estimation.
# Updated: 2026-06-16

#' Calculate AUC Components for One Arm and Event Type
#'
#' @param data Per-process data for a single arm.
#' @param tau Truncation time.
#' @param calc_var Calculate variance components?
#' @return List with area, mcf, psi, n, var_area.
#' @noRd
MvProcessAUC <- function(
    data,
    tau,
    calc_var = TRUE
) {

  if (NROW(data) < 1) {
    stop("No observations available for the requested arm and event type.")
  }

  n <- length(unique(data$idx))
  if (!any(data$status == 1)) {
    out <- list(
      area = 0,
      mcf = data.frame(
        time = 0,
        censor = 0,
        death = 0,
        events = 0,
        nar = n,
        haz = 0,
        surv = 1,
        event_rate = 0,
        weighted_event_rate = 0,
        mcf = 0,
        var_mcf = 0,
        se_mcf = 0
      ),
      n = n
    )
    if (calc_var) {
      out$psi <- data.frame(
        idx = unique(data$idx),
        psi = 0
      )
      out$var_area <- 0
    }
    return(out)
  }

  mcf <- CalcMCF(
    idx = data$idx,
    status = data$status,
    time = data$time,
    jump_weights = data$jump_weights,
    calc_var = calc_var
  )

  area <- AUC(
    times = mcf$time,
    values = mcf$mcf,
    tau = tau
  )

  n <- length(unique(data$idx))
  out <- list(
    area = area,
    mcf = mcf,
    n = n
  )

  if (calc_var) {
    psi_df <- VarAUC(
      data = data[, c("idx", "status", "time", "jump_weights")],
      tau = tau,
      mcf = mcf,
      return_psi = TRUE,
      jump_weights = data$jump_weights
    )
    out$psi <- psi_df
    out$var_area <- mean(psi_df$psi^2)
  }

  return(out)
}


#' Calculate Multivariate Theta and Influence Functions for One Arm
#'
#' @param arm_data Arm-specific formatted data.
#' @param event_types Event types.
#' @param tau Truncation time.
#' @param cens_after_last Censor after last event?
#' @param covar_cols Covariate column names, or NULL.
#' @return List with per-process results and stacked influence matrix.
#' @noRd
MvArmTheta <- function(
    arm_data,
    event_types,
    tau,
    cens_after_last = TRUE,
    covar_cols = NULL
) {

  K <- length(event_types)
  arm <- unique(arm_data$arm)
  all_idx <- sort(unique(arm_data$idx))
  n_arm <- length(all_idx)

  areas <- numeric(K)
  names(areas) <- as.character(event_types)
  mcf_list <- vector("list", K)
  names(mcf_list) <- as.character(event_types)
  n_ak <- numeric(K)
  names(n_ak) <- as.character(event_types)
  psi_mat <- matrix(0, nrow = n_arm, ncol = K)
  rownames(psi_mat) <- as.character(all_idx)
  colnames(psi_mat) <- as.character(event_types)
  xbar_by_k <- vector("list", K)

  for (k in seq_len(K)) {
    et <- event_types[k]
    at_risk_idx <- .MvAtRiskIdx(data = arm_data, k = et)

    if (length(at_risk_idx) == 0) {
      stop(
        "No at-risk subjects for arm ", arm, " and event type ", et, "."
      )
    }

    arm_k <- arm_data[arm_data$idx %in% at_risk_idx, , drop = FALSE]
    proc <- .SubsetProcessK(
      data = arm_k,
      k = et,
      cens_after_last = cens_after_last
    )

    fit <- MvProcessAUC(
      data = proc,
      tau = tau,
      calc_var = TRUE
    )
    areas[k] <- fit$area
    mcf_list[[as.character(et)]] <- fit$mcf
    n_ak[k] <- fit$n

    psi_k <- fit$psi
    eligibility_weight <- n_arm / fit$n

    row_match <- match(as.character(psi_k$idx), rownames(psi_mat))
    psi_mat[row_match, k] <- eligibility_weight * psi_k$psi

    if (!is.null(covar_cols)) {
      cov_sub <- arm_k %>%
        dplyr::distinct(idx, dplyr::across(dplyr::all_of(covar_cols))) %>%
        dplyr::filter(idx %in% at_risk_idx)
      xbar_by_k[[k]] <- colMeans(
        as.matrix(cov_sub[, covar_cols, drop = FALSE]),
        na.rm = TRUE
      )
    }
  }

  out <- list(
    arm = arm,
    areas = areas,
    mcf_list = mcf_list,
    n_arm = n_arm,
    n_ak = n_ak,
    psi_mat = psi_mat,
    idx = all_idx,
    xbar_by_k = xbar_by_k
  )
  return(out)
}


#' Calculate Multivariate Augmented AUC Statistics
#'
#' @param data Formatted multivariate data.
#' @param event_types Event types.
#' @param tau Truncation time.
#' @param alpha Type I error level.
#' @param cens_after_last Censor after last event?
#' @param covar_cols Covariate column names or NULL.
#' @param process_weights Optional contrast weights of length K.
#' @param return_areas Return area summaries?
#' @return List of test statistics and covariance components.
#' @noRd
CalcMvAugAUC <- function(
    data,
    event_types,
    tau,
    alpha = 0.05,
    cens_after_last = TRUE,
    covar_cols = NULL,
    process_weights = NULL,
    return_areas = FALSE
) {

  K <- length(event_types)

  data0 <- data[data[["arm"]] == 0, , drop = FALSE]
  data1 <- data[data[["arm"]] == 1, , drop = FALSE]

  a0 <- MvArmTheta(
    arm_data = data0,
    event_types = event_types,
    tau = tau,
    cens_after_last = cens_after_last,
    covar_cols = covar_cols
  )
  a1 <- MvArmTheta(
    arm_data = data1,
    event_types = event_types,
    tau = tau,
    cens_after_last = cens_after_last,
    covar_cols = covar_cols
  )

  theta0 <- a0$areas
  theta1 <- a1$areas
  delta <- theta1 - theta0

  sigma0 <- crossprod(a0$psi_mat) / a0$n_arm
  sigma1 <- crossprod(a1$psi_mat) / a1$n_arm
  cov_delta <- sigma1 / a1$n_arm + sigma0 / a0$n_arm

  if (!is.null(covar_cols)) {
    cov0 <- data0 %>%
      dplyr::distinct(idx, dplyr::across(dplyr::all_of(covar_cols))) %>%
      dplyr::arrange(idx)
    cov1 <- data1 %>%
      dplyr::distinct(idx, dplyr::across(dplyr::all_of(covar_cols))) %>%
      dplyr::arrange(idx)
    cov0 <- as.matrix(cov0[match(a0$idx, cov0$idx), covar_cols, drop = FALSE])
    cov1 <- as.matrix(cov1[match(a1$idx, cov1$idx), covar_cols, drop = FALSE])

    aug0 <- CalcMvAugCompArm(
      covars = cov0,
      psi_mat = a0$psi_mat,
      xbar_by_k = a0$xbar_by_k,
      n_a = a0$n_arm
    )
    aug1 <- CalcMvAugCompArm(
      covars = cov1,
      psi_mat = a1$psi_mat,
      xbar_by_k = a1$xbar_by_k,
      n_a = a1$n_arm
    )

    sigma_imb <- aug0$sigma + aug1$sigma
    gamma <- aug0$gamma + aug1$gamma

    p <- ncol(cov0)
    delta_z <- numeric(K * p)
    for (k in seq_len(K)) {
      cols <- ((k - 1L) * p + 1L):(k * p)
      delta_z[cols] <- a1$xbar_by_k[[k]] - a0$xbar_by_k[[k]]
    }

    B <- MASS::ginv(sigma_imb) %*% gamma
    correction <- as.numeric(crossprod(B, delta_z))
    delta <- delta - correction
    cov_delta <- cov_delta - crossprod(B, gamma)
  }

  se_delta <- sqrt(pmax(diag(cov_delta), 0))
  crit <- stats::qnorm(1 - alpha / 2)

  contrasts <- data.frame(
    contrast = "A1-A0",
    event_type = event_types,
    observed = as.numeric(delta),
    se = se_delta,
    lower = as.numeric(delta) - crit * se_delta,
    upper = as.numeric(delta) + crit * se_delta,
    p = 2 * stats::pnorm(abs(as.numeric(delta)) / se_delta, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )

  weighted <- data.frame()
  if (!is.null(process_weights)) {
    if (length(process_weights) != K) {
      stop("process_weights must have length K = ", K, ".")
    }
    w <- as.numeric(process_weights)
    obs_w <- sum(w * delta)
    var_w <- as.numeric(t(w) %*% cov_delta %*% w)
    se_w <- sqrt(max(var_w, 0))
    weighted <- data.frame(
      contrast = "w'Delta",
      observed = obs_w,
      se = se_w,
      lower = obs_w - crit * se_w,
      upper = obs_w + crit * se_w,
      p = 2 * stats::pnorm(abs(obs_w) / se_w, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    contrasts = contrasts,
    cov_delta = cov_delta,
    cov_theta = list("0" = sigma0, "1" = sigma1),
    theta0 = theta0,
    theta1 = theta1,
    weighted = weighted,
    a0 = a0,
    a1 = a1,
    tau = tau
  )

  if (return_areas) {
    se0 <- sqrt(diag(sigma0) / a0$n_arm)
    se1 <- sqrt(diag(sigma1) / a1$n_arm)
    areas0 <- data.frame(
      arm = 0,
      event_type = event_types,
      n = a0$n_ak,
      tau = tau,
      area = as.numeric(theta0),
      se = se0
    )
    areas1 <- data.frame(
      arm = 1,
      event_type = event_types,
      n = a1$n_ak,
      tau = tau,
      area = as.numeric(theta1),
      se = se1
    )
    out$marg_areas <- rbind(areas0, areas1)

    mcf0 <- lapply(names(a0$mcf_list), function(et) {
      m <- a0$mcf_list[[et]]
      m$arm <- 0
      m$event_type <- as.integer(et)
      return(m)
    })
    mcf1 <- lapply(names(a1$mcf_list), function(et) {
      m <- a1$mcf_list[[et]]
      m$arm <- 1
      m$event_type <- as.integer(et)
      return(m)
    })
    out$mcf <- do.call(rbind, c(mcf0, mcf1))
  }

  return(out)
}
