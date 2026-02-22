## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache = TRUE)
library(dplyr)
library(ggplot2)
library(MCC)

## -----------------------------------------------------------------------------
N <- 1e3
DGP <- function() {
  out <- GenData(
    base_event_rate = 0.5,
    base_death_rate = 0.2,
    censoring_rate = 0.2,
    tau = 5,
    n = N
  )
  out$weights <- 2
  return(out)
}

## -----------------------------------------------------------------------------
#' Calculate MCF
GetMCF <- function(df) {
  mcf <- CalcMCF(
    idx = df$idx,
    status = df$status,
    time = df$time,
    weights = df$weights
  )
  return(mcf)
}

#' MCF curve.
Curve <- function(mcf) {
  out <- stats::stepfun(x = mcf$time, y = c(0, mcf$mcf))
  return(out)
}

#' SE curve.
SeCurve <- function(mcf) {
  out <- stats::stepfun(x = mcf$time, y = c(0, mcf$se_mcf))
  return(out)
}

#' Simulation loop.
Loop <- function(i) {
  
  # Generate data.
  data <- DGP()
  
  # Calculate MCF.
  mcf <- GetMCF(data)
  mcf_fn <- Curve(mcf)
  mcf_se <- SeCurve(mcf)  
  
  # Evaluate.
  taus <- seq(1:4)
  mcf_evals <- sapply(taus, mcf_fn)
  se_evals <- sapply(taus, mcf_se)
  
  # Output.
  out <- data.frame(
    idx = i,
    tau = taus,
    mcf = mcf_evals,
    ses = se_evals
  )
  return(out)
}

#' Simulation.
Sim <- function(reps = 1e3) {
  
  sim <- lapply(seq_len(reps), Loop)
  sim <- do.call(rbind, sim)
  
  results <- sim %>%
    dplyr::group_by(tau) %>%
    dplyr::summarise(
      reps = dplyr::n(),
      empirical_se = stats::sd(mcf),
      analytical_se = sqrt(mean(ses^2))
    )
  
  return(results)
}


## ----cache=TRUE---------------------------------------------------------------
sim <- Sim(reps = 500)

## ----echo = FALSE, fig.align='center', cache=TRUE-----------------------------
plot_df <- sim %>%
  dplyr::select(tau, empirical_se, analytical_se) %>%
  tidyr::pivot_longer(
    empirical_se:analytical_se,
    names_to = "method",
    values_to = "se"
  )

plot_df$method <- factor(
  x = plot_df$method, 
  levels = c("empirical_se", "analytical_se"),
  labels = c("Empirical", "Analytical")
)

plot_df$scaled_se <- plot_df$se * sqrt(N)

q <- ggplot(data = plot_df) + 
  theme_bw() + 
  theme(legend.position = "top") + 
  geom_col(
    aes(x = tau, y = scaled_se, fill = method),
    position = position_dodge()
  ) + 
  ggsci::scale_fill_nejm(name = "Method") + 
  scale_x_continuous(name = "Truncation Time") + 
  scale_y_continuous(name = expression(sqrt(N) %*% "SE"), limits = c(0, 4.0))
show(q)

## -----------------------------------------------------------------------------
N <- 250
DGP <- function() {
  out <- GenData(
    base_event_rate = 0.5,
    base_death_rate = 0.2,
    censoring_rate = 0.2,
    tau = 5,
    n = N
  )
  out$weights <- 2
  return(out)
}

## -----------------------------------------------------------------------------
#' Calculate MCF
GetAUC <- function(df, tau) {
  auc <- SingleArmAUC(
    data = df,
    weights = df$weights,
    tau = tau
  )
  out <- auc@MargAreas
  return(out)
}

#' Simulation loop.
Loop <- function(i) {
  
  # Generate data.
  data <- DGP()
  
  # Calculate MCF.
  taus <- c(1, 2, 3, 4)
  out <- lapply(taus, function(t) {GetAUC(data, t)})
  out <- data.frame(do.call(rbind, out))
  
  # Output.
  out$arm <- NULL
  out$idx <- i
  return(out)
}

#' Simulation.
Sim <- function(reps = 1e3) {
  
  sim <- lapply(seq_len(reps), Loop)
  sim <- do.call(rbind, sim)
  
  results <- sim %>%
    dplyr::group_by(tau) %>%
    dplyr::summarise(
      reps = dplyr::n(),
      empirical_se = stats::sd(area),
      analytical_se = sqrt(mean(se^2))
    )
  
  return(results)
}


## ----cache=TRUE---------------------------------------------------------------
sim <- Sim(reps = 500)

## ----echo = FALSE, fig.align='center', cache=TRUE-----------------------------
plot_df <- sim %>%
  dplyr::select(tau, empirical_se, analytical_se) %>%
  tidyr::pivot_longer(
    empirical_se:analytical_se,
    names_to = "method",
    values_to = "se"
  )

plot_df$method <- factor(
  x = plot_df$method, 
  levels = c("empirical_se", "analytical_se"),
  labels = c("Empirical", "Analytical")
)

plot_df$scaled_se <- plot_df$se * sqrt(N)

q <- ggplot(data = plot_df) + 
  theme_bw() + 
  theme(legend.position = "top") + 
  geom_col(
    aes(x = tau, y = scaled_se, fill = method),
    position = position_dodge()
  ) + 
  ggsci::scale_fill_nejm(name = "Method") + 
  scale_x_continuous(name = "Truncation Time") + 
  scale_y_continuous(name = expression(sqrt(N) %*% "SE"), limits = c(0, 8.0))
show(q)

