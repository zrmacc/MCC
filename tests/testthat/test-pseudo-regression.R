.reference_h12 <- function(prepared, design, tau, integrate = FALSE) {
  d <- prepared$process_data
  grid <- prepared$mcf$time[prepared$mcf$time <= tau]
  ids <- sort(unique(d$idx))
  n <- length(ids)
  m <- length(grid)
  risk <- events <- deaths <- matrix(0, m, n)
  for (i in seq_along(ids)) {
    di <- d[d$idx == ids[i], , drop = FALSE]
    risk[, i] <- as.numeric(grid <= max(di$time))
    for (r in seq_len(nrow(di))) {
      u <- match(di$time[r], grid)
      if (is.na(u)) next
      if (di$status[r] == 1) events[u, i] <- events[u, i] + di$jump_weights[r]
      if (di$status[r] == 2) deaths[u, i] <- 1
    }
  }
  y <- rowMeans(risk)
  da <- rowMeans(events)
  da_d <- rowMeans(deaths)
  dl <- da_d / y
  denom <- 1 - dl
  s_left <- c(1, utils::head(prepared$mcf$surv, m - 1L))
  b <- s_left / y
  z <- lapply(seq_len(n), function(i) list(y = risk[, i] - y, a = events[, i] - da, d = deaths[, i] - da_d))
  first <- function(h) {
    dl_h <- h$d / y - h$y * da_d / y^2
    j_left <- c(0, utils::head(cumsum(dl_h / denom), -1L))
    list(dl = dl_h, r = j_left + h$y / y)
  }
  mixed <- function(h, k) {
    fh <- first(h)
    fk <- first(k)
    dl_hk <- -h$y * k$d / y^2 - k$y * h$d / y^2 + 2 * h$y * k$y * da_d / y^3
    j_hk_left <- c(0, utils::head(cumsum(dl_hk / denom + fh$dl * fk$dl / denom^2), -1L))
    r_hk <- j_hk_left - h$y * k$y / y^2
    increment <- b * (fh$r * fk$r - r_hk) * da - b * fh$r * k$a - b * fk$r * h$a
    sum(increment * if (integrate) pmax(tau - grid, 0) else 1)
  }
  out <- matrix(0, n, ncol(design), dimnames = list(NULL, colnames(design)))
  for (i in seq_len(n)) for (k in seq_len(ncol(design))) {
    out[i, k] <- mean(vapply(seq_len(n), function(j) design[j, k] * mixed(z[[j]], z[[i]]), numeric(1)))
  }
  out
}

test_that("PseudoReg matches lm coefficients and manual corrected covariance", {
  withr::local_seed(901)
  x <- rnorm(30)
  w <- cbind(intercept = 1, A = rep(0:1, each = 15), X = x)
  data <- GenData(covariates = w, beta_event = c(0, .2, -.1), beta_death = c(0, 0, 0), n = 30, tau = 3)
  fit <- PseudoReg(~ A + X + A:X, data, tau = 2, type = "MCF")
  plain <- lm(pseudo ~ A + X + A:X, data = cbind(fit$pseudo, A = w[, "A"], X = w[, "X"]))
  expect_equal(coef(fit), coef(plain), tolerance = 1e-12)
  expect_lt(max(abs(colMeans(fit$h12))), 1e-12)
  centered <- sweep(fit$score, 2, colMeans(fit$score), "-")
  manual <- solve(fit$Ahat) %*% (crossprod(centered) / nrow(centered)) %*% solve(fit$Ahat) / nrow(centered)
  expect_equal(vcov(fit), manual, tolerance = 1e-12)
  expect_equal(confint(fit), cbind(fit$coefficients$lower, fit$coefficients$upper), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("optimized h12 matches the direct pairwise reference", {
  withr::local_seed(902)
  w <- cbind(intercept = 1, A = rep(0:1, each = 6), X = rnorm(12))
  data <- GenData(covariates = w, beta_event = c(0, .2, .1), beta_death = c(0, 0, 0), n = 12, tau = 3)
  for (type in c("MCF", "AUC")) {
    fit <- PseudoReg(~ A + X, data, tau = 2, type = type)
    prepared <- MCC:::.PreparePseudo(data, tau = 2, type = type)
    reference <- .reference_h12(prepared, fit$fit$x, tau = 2, integrate = type == "AUC")
    expect_equal(unname(fit$h12), unname(reference), tolerance = 2e-11)
  }
})

test_that("intercept-only and no-censoring projections vanish", {
  withr::local_seed(903)
  w <- cbind(intercept = 1, X = rnorm(20))
  censored <- GenData(covariates = w, beta_event = c(0, .1), beta_death = c(0, 0), n = 20, tau = 3)
  expect_lt(max(abs(PseudoReg(~ 1, censored, tau = 2)$h12)), 1e-12)
  data <- GenData(covariates = w, beta_event = c(0, .1), beta_death = c(0, 0), censoring_rate = 0, n = 20, tau = 3)
  fit <- PseudoReg(~ X, data, tau = 2)
  expect_lt(max(abs(fit$h12)), 1e-12)
  bread <- solve(crossprod(fit$fit$x))
  hc0 <- bread %*% crossprod(fit$fit$x * residuals(fit$fit)) %*% bread
  expect_equal(vcov(fit), hc0, tolerance = 1e-12)
})

test_that("empirical derivatives agree with centered finite differences", {
  withr::local_seed(906)
  w <- cbind(intercept = 1, X = rnorm(10))
  data <- GenData(covariates = w, beta_event = c(0, .1), beta_death = c(0, 0), n = 10, tau = 3)
  prepared <- MCC:::.PreparePseudo(data, tau = 2, type = "MCF")
  d <- prepared$process_data
  grid <- prepared$mcf$time[prepared$mcf$time <= 2]
  ids <- sort(unique(d$idx))
  n <- length(ids)
  m <- length(grid)
  risk <- events <- deaths <- matrix(0, m, n)
  for (i in seq_along(ids)) {
    di <- d[d$idx == ids[i], , drop = FALSE]
    risk[, i] <- as.numeric(grid <= max(di$time))
    for (r in seq_len(nrow(di))) {
      u <- match(di$time[r], grid)
      if (is.na(u)) next
      if (di$status[r] == 1) events[u, i] <- events[u, i] + di$jump_weights[r]
      if (di$status[r] == 2) deaths[u, i] <- 1
    }
  }
  base <- list(y = rowMeans(risk), a = rowMeans(events), d = rowMeans(deaths))
  z <- lapply(seq_len(n), function(i) list(y = risk[, i] - base$y, a = events[, i] - base$a, d = deaths[, i] - base$d))
  functional <- function(e1, e2, h, k) {
    y <- base$y + e1 * h$y + e2 * k$y
    da <- base$a + e1 * h$a + e2 * k$a
    da_d <- base$d + e1 * h$d + e2 * k$d
    survival <- cumprod(1 - da_d / y)
    survival_left <- c(1, utils::head(survival, -1L))
    sum(survival_left / y * da)
  }
  eps1 <- 1e-6
  first_fd <- (functional(eps1, 0, z[[2]], z[[3]]) - functional(-eps1, 0, z[[2]], z[[3]])) / (2 * eps1)
  dl <- base$d / base$y
  dl_h <- z[[2]]$d / base$y - z[[2]]$y * base$d / base$y^2
  j_h_left <- c(0, utils::head(cumsum(dl_h / (1 - dl)), -1L))
  r_h <- j_h_left + z[[2]]$y / base$y
  survival_left <- c(1, utils::head(cumprod(1 - dl), -1L))
  b <- survival_left / base$y
  first_analytic <- sum(-b * r_h * base$a + b * z[[2]]$a)
  expect_equal(first_fd, first_analytic, tolerance = 2e-7)
  eps2 <- 2e-4
  mixed_fd <- (functional(eps2, eps2, z[[2]], z[[3]]) - functional(eps2, -eps2, z[[2]], z[[3]]) - functional(-eps2, eps2, z[[2]], z[[3]]) + functional(-eps2, -eps2, z[[2]], z[[3]])) / (4 * eps2^2)
  pair_design <- diag(n)
  direct_all <- .reference_h12(prepared, pair_design, tau = 2)
  expect_equal(mixed_fd, n * direct_all[3, 2], tolerance = 2e-6)
})

test_that("PseudoReg supports separate factors and validates subject alignment", {
  withr::local_seed(904)
  w <- cbind(intercept = 1, A = rep(0:1, each = 10), X = rnorm(20))
  data <- GenData(covariates = w, beta_event = c(0, 0, 0), beta_death = c(0, 0, 0), n = 20, tau = 3)
  subject <- data.frame(idx = seq_len(20), group = factor(rep(c("a", "b"), each = 10)), X = w[, "X"])
  fit <- PseudoReg(~ group * X, data, covariates = subject, tau = 2)
  expect_s3_class(fit, "PseudoReg")
  expect_equal(colnames(fit$h12), names(coef(fit)))
  expect_error(PseudoReg(~ group, data, covariates = subject[-1, ], tau = 2), "do not match")
  duplicate <- rbind(subject, subject[1, ])
  expect_error(PseudoReg(~ group, data, covariates = duplicate, tau = 2), "one row per subject")
  varying <- data
  varying$X[varying$idx == 1][1] <- varying$X[varying$idx == 1][1] + 1
  expect_error(PseudoReg(~ X, varying, tau = 2), "not constant")
})

test_that("PseudoReg rejects unsupported model and data configurations", {
  withr::local_seed(905)
  w <- cbind(intercept = 1, X = rnorm(12))
  data <- GenData(covariates = w, beta_event = c(0, 0), beta_death = c(0, 0), n = 12, tau = 3)
  expect_error(PseudoReg(pseudo ~ X, data, tau = 2), "one-sided")
  expect_error(PseudoReg(~ ., data, tau = 2), "shorthand")
  expect_error(PseudoReg(~ offset(X), data, tau = 2), "Offsets")
  data$copy <- data$X
  expect_error(PseudoReg(~ X + copy, data, tau = 2), "rank deficient")
  expect_error(PseudoReg(~ X, data, tau = 4), "observed time range")
  shifted <- data
  shifted$time <- shifted$time + 1
  expect_error(PseudoReg(~ X, shifted, tau = 0.5), "observed time range")
  exhausted <- data.frame(idx = rep(1:4, each = 2), time = rep(c(.5, 1), 4), status = rep(c(1, 2), 4), X = rep(1:4, each = 2))
  expect_error(PseudoReg(~ X, exhausted, tau = 1), "risk set is exhausted")
})

test_that("custom column names and jump weights are supported", {
  withr::local_seed(907)
  w <- cbind(intercept = 1, X = rnorm(16))
  data <- GenData(covariates = w, beta_event = c(0, .1), beta_death = c(0, 0), n = 16, tau = 3)
  names(data)[match(c("idx", "status", "time"), names(data))] <- c("patient", "code", "when")
  weights <- ifelse(data$code == 1, 1.5, 1)
  fit <- PseudoReg(~ X, data, tau = 1.75, idx_name = "patient", status_name = "code", time_name = "when", jump_weights = weights)
  expect_true(all(is.finite(fit$vcov)))
  expect_equal(fit$tau, 1.75)
})

test_that("nonconstant jump weights remain aligned when input rows are reordered", {
  withr::local_seed(908)
  n <- 20
  x <- rnorm(n)
  w <- cbind(intercept = 1, A = rep(0:1, each = n / 2), X = x)
  data <- GenData(covariates = w, beta_event = c(0, .2, -.1), beta_death = c(0, 0, 0), n = n, tau = 3)
  subject <- data.frame(idx = seq_len(n), A = w[, "A"], X = x)
  weights <- ifelse(data$status == 1, seq_len(nrow(data)) / nrow(data) + .5, 1)
  original <- PseudoReg(~ A + X, data, covariates = subject, tau = 2, jump_weights = weights)
  row_order <- sample.int(nrow(data))
  subject_order <- sample.int(n)
  reordered <- PseudoReg(
    ~ A + X,
    data[row_order, , drop = FALSE],
    covariates = subject[subject_order, , drop = FALSE],
    tau = 2,
    jump_weights = weights[row_order]
  )
  expect_equal(coef(reordered), coef(original), tolerance = 1e-12)
  expect_equal(vcov(reordered), vcov(original), tolerance = 1e-12)
  expect_equal(
    reordered$h12[order(rownames(reordered$h12)), , drop = FALSE],
    original$h12[order(rownames(original$h12)), , drop = FALSE],
    tolerance = 1e-12
  )
})
