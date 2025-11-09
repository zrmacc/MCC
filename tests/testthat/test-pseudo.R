test_that("Compare pseudo-value calculation with jackknife.", {
  
  AUC <- function(data, tau) {
    fit <- SingleArmAUC(data = data, tau = tau)
    return(fit@MargAreas$area)
  }
  
  Jackknife <- function(data, tau) {
    orig <- AUC(data, tau)
    n <- length(unique(data$idx))
    out <- lapply(seq_len(n), function(i) {
      sub <- data[data$idx != i, ]
      sub_auc <- AUC(data = sub, tau = tau)
      out <- data.frame(
        idx = i,
        jk = sub_auc + n * (orig - sub_auc)
      )
      return(out)
    })
    out <- do.call(rbind, out)
    return(out)
  }
  
  data <- data.frame(
    idx = c(1, 1, 2, 2, 3, 3, 4, 4),
    time = c(1.0, 2.5, 1.2, 3.0, 0.5, 2.0, 1.8, 4.0),
    status = c(1, 0, 1, 0, 1, 0, 1, 0) 
  )
  
  # tau = 1.0
  taus <- seq(1:4)
  for (tau in taus) {
    obs <- GenPseudo(data = data, tau = tau)
    exp <- Jackknife(data = data, tau = tau)  
    expect_equal(obs$pseudo, exp$jk)
  }
  
})


# -----------------------------------------------------------------------------

test_that("Overall test of pseudo-values.", {
  
  # Calculate difference of AUCs via package.
  DiffAUCs <- function(data, tau) {
    comp <- CompareAUCs(data = data, tau = tau)
    out <- comp@Pvals$observed[1]
    return(out)
  }
  
  # Calculate difference of AUCs via linear regression.
  LM <- function(data, tau) {
    df_pseudo <- GenPseudo(data = data, tau = tau)
    df_pseudo$arm <- covar$arm[df_pseudo$idx]
    fit <- lm(pseudo ~ arm, df_pseudo)
    out <- as.numeric(coef(fit)[2])
    return(out)
  }
  
  set.seed(101)
  covar <- data.frame(arm = c(rep(1, 50), rep(0, 50)))
  data <- GenData(beta_event = log(0.5), tau = 4, covar = covar)

  taus <- seq(1:4)
  for (tau in taus) {
    obs <- DiffAUCs(data = data, tau = tau)
    exp <- LM(data = data, tau = tau)  
    expect_equal(obs, exp, tolerance = 0.01)
  }
  
})

