
# Comparison of mean cumulative count curves via the area under the curve

[![R-CMD-check](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml)

Zachary R. McCaw <br> Updated: 2026-06-25

## Installation

``` r
remotes::install_github("zrmacc/MCC", build_vignettes = TRUE)
```

## Description

This package provides estimation and inference for the mean cumulative
function (MCF) and the area under the MCF (AUMCF). The mean cumulative
count (MCC) curve is an alias for the MCF. Estimation and inference for
the MCF follow [Ghosh and Lin
(2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x);
corresponding procedures for the AUMCF are described by [Gronsbell et
al. (2026)](https://academic.oup.com/biometrics/article/82/2/ujag053/8689483).

## Univariate

Univariate analysis considers a single multiple- or recurrent-events
process, subject to right-censoring and terminal events.

### Methods

#### MCF estimation

For each study arm, define $N(\cdot)$ as the counting process for events
of interest (terminal and non-terminal), $Y(\cdot)$ as the number of
subjects at risk, and $S(\cdot)$ as the probability of not having
experienced a terminal event. The MCF $\mu(t)$ at time $t$ is estimated
by:

$$
\mu(t) = \int_{0}^{t}\ \hat{S}(u)\ \frac{ dN(u) }{ Y(u) }
$$

Here $\hat{S}(u)$ is the Kaplan-Meier estimate of the probability of
being terminal event-free; $dN(u)$ is the number of events of interest
at time $u$; and $Y(u)$ is the number of subjects who remain at risk.
The AUMCF to time $\tau$ is
$\hat\alpha(\tau) = \int_0^\tau \hat\mu(u)\,du$.

#### Influence function

Let $\theta(\tau)$ denote the MCF or AUMCF at time $\tau$. The influence
function $\psi_{i}(\tau)$ measures the contribution of subject $i$ to
the estimation error:

$$
\hat{\theta}(\tau) - \theta(\tau) = \frac{1}{n}\sum_{i=1}^{n}\psi_{i}(\tau) + o_{p}(n^{-1/2})
$$

The influence function contributions estimate the variance of
$\hat{\theta}(\tau)$:

$$
\hat{\mathbb{V}}(\hat{\theta}) = \frac{1}{n^2}\sum_{i=1}^{n}\psi_{i}^{2}(\tau).
$$

### Data

The function `GenData` simulates example data. Recurrent event times are
generated from a Poisson process until censoring or death. Optionally, a
gamma `frailty_variance` correlates patient-specific event and death
rates. The example below has 100 patients per arm, `tau = 4`, and an 80%
recurrent-event rate in the treatment arm relative to control.

``` r
covariates <- data.frame(
  arm = c(rep(1, 100), rep(0, 100))
)
data <- MCC::GenData(
  beta_event = c(log(0.8)),
  covariates = covariates,
  frailty_variance = 0.2,
  tau = 4
)
head(data)
```

    ##   idx status      time arm cens_rate death_rate event_rate   frailty
    ## 1   1      0 0.5900957   1      0.25  0.3946803  1.2629769 1.5787211
    ## 2   2      1 0.2753939   1      0.25  0.1688872  0.5404391 0.6755488
    ## 3   2      1 0.6651751   1      0.25  0.1688872  0.5404391 0.6755488
    ## 4   2      0 3.0821456   1      0.25  0.1688872  0.5404391 0.6755488
    ## 5   3      1 2.3232824   1      0.25  0.1282805  0.4104975 0.5131219
    ## 6   3      0 4.0000000   1      0.25  0.1282805  0.4104975 0.5131219

The essential columns are:

- `idx`, the subject index.
- `time`, the observation time.
- `status`, coded 0 for censoring, 1 for an event, 2 for death (or any
  competing terminal event).
- `arm`, coded as 1 for treatment, 0 for reference.

Each subject should experience an observation-terminating event
(censoring or death). For other data sets, use the same arm and status
coding.

#### Observation-terminating events

In the recurrent-events setting, a subject may remain at risk after an
event of interest. An *observation-terminating* event (censoring or
competing risk) removes the subject from the risk set. If a subject
lacks such an event, they are assumed to remain at risk indefinitely. By
default, `CompareAUCs` adds a censoring time immediately after the last
event of interest. For example, if the data for subject 1 were:

    ##   idx time status
    ## 1   1    2      1
    ## 2   1    3      1
    ## 3   1    5      1

then, for analysis, the subject is assumed censored after the last
event:

    ##   idx time status
    ## 1   1    2      1
    ## 2   1    3      1
    ## 3   1    5      1
    ## 4   1    5      0

Set `cens_after_last = FALSE` if the subject should remain at risk
indefinitely.

#### Terminal events of interest

When a fatal event is of interest, record it with two rows: `status = 1`
(event occurred) and `status = 2` (event was terminal). For example,
subject 1 had three events and the third was terminal:

    ##   idx time status
    ## 1   1    2      1
    ## 2   1    3      1
    ## 3   1    5      1
    ## 4   1    5      2

By contrast, subject 2 had three non-terminal events and is censored
after the last by default:

    ##   idx time status
    ## 1   2    2      1
    ## 2   2    3      1
    ## 3   2    5      1
    ## 4   2    5      0

Censoring (`status = 0`) and terminal events (`status = 2`) both remove
a subject from the risk set, but only censoring leaves open the
possibility of future events.

### Analyses

#### Single-arm AUC

``` r
auc <- MCC::SingleArmAUC(
  data %>% dplyr::filter(arm == 0),
  boot = TRUE,
  reps = 200,
  tau = 4
)
show(auc)
```

    ## Marginal Areas:
    ##   arm   n area    se tau
    ## 1   0 100 6.65 0.633   4
    ## 
    ## 
    ## CIs:
    ##       method contrast observed    se lower upper
    ## 1 asymptotic       A0     6.65 0.633  5.41  7.89
    ## 2  bootstrap       A0     6.65 0.632  5.32  7.72
    ## 
    ## 
    ## P-values:
    ##       method contrast observed        p
    ## 1 asymptotic       A0     6.65 8.21e-26

#### Two-arm comparison

``` r
aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  boot = TRUE,
  perm = TRUE,
  reps = 200,
  alpha = 0.05
)
show(aucs)
```

    ## Marginal Areas:
    ##   arm   n area    se tau
    ## 1   0 100 6.65 0.633   4
    ## 2   1 100 4.26 0.516   4
    ## 
    ## 
    ## CIs:
    ##       method contrast observed     se  lower  upper
    ## 1 asymptotic    A1-A0    -2.39 0.8160 -3.990 -0.789
    ## 3  bootstrap    A1-A0    -2.39 0.7630 -3.600 -0.812
    ## 2 asymptotic    A1/A0     0.64 0.0987  0.473  0.866
    ## 4  bootstrap    A1/A0     0.64 0.0971  0.487  0.854
    ## 
    ## 
    ## P-values:
    ##        method contrast observed       p
    ## 1  asymptotic    A1-A0    -2.39 0.00343
    ## 3   bootstrap    A1-A0    -2.39 0.01990
    ## 5 permutation    A1-A0    -2.39 0.02990
    ## 2  asymptotic    A1/A0     0.64 0.00385
    ## 4   bootstrap    A1/A0     0.64 0.01990
    ## 6 permutation    A1/A0     0.64 0.02990

- `tau` is the truncation time for the AUC.
- `boot` constructs bootstrap confidence intervals.
- `perm` performs permutation tests for the difference and ratio of
  AUCs.
- `reps` is the number of bootstrap/permutation replicates (grouped by
  `idx`, stratified by `strata` if applicable).
- `alpha` is 1 minus the desired CI coverage.

#### Jump-weighted analysis

The `jump_weights` argument scales the jump in the MCF at each recurrent
event (`status == 1`). The example below weights each event by the
patient’s cumulative event count:

``` r
data <- data %>%
  dplyr::group_by(idx) %>%
  dplyr::mutate(jump_weights = dplyr::row_number()) %>%
  dplyr::ungroup()

cat("Jump weights for the first 10 records.\n")
data %>%
  dplyr::select(idx, time, status, jump_weights) %>%
  dplyr::slice(1:10)
```

    ## Jump weights for the first 10 records.
    ## # A tibble: 10 × 4
    ##      idx  time status jump_weights
    ##    <dbl> <dbl>  <dbl>        <int>
    ##  1     1 0.341      1            1
    ##  2     1 1.67       0            2
    ##  3     2 2.54       1            1
    ##  4     2 3.16       1            2
    ##  5     2 4          0            3
    ##  6     3 0.434      1            1
    ##  7     3 0.571      1            2
    ##  8     3 1.49       1            3
    ##  9     3 1.86       0            4
    ## 10     4 0.597      2            1

``` r
aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  alpha = 0.05,
  jump_weights = data$jump_weights
)
show(aucs)
```

    ## Marginal Areas:
    ##   arm   n area   se tau
    ## 1   0 100 16.5 2.62   4
    ## 2   1 100  8.6 1.51   4
    ## 
    ## 
    ## CIs:
    ##       method contrast observed    se   lower  upper
    ## 1 asymptotic    A1-A0   -7.860 3.030 -13.800 -1.930
    ## 2 asymptotic    A1/A0    0.522 0.124   0.328  0.831
    ## 
    ## 
    ## P-values:
    ##       method contrast observed       p
    ## 1 asymptotic    A1-A0   -7.860 0.00934
    ## 2 asymptotic    A1/A0    0.522 0.00614

#### Stratified analysis

The example below adds a binary stratification factor; the event rate in
stratum 1 is increased by 20%.

``` r
covariates <- data.frame(
  arm = c(rep(1, 100), rep(0, 100)),
  strata = stats::rbinom(200, 1, 0.25)
)
data <- MCC::GenData(
  beta_event = c(log(0.8), log(1.2)),
  covariates = covariates,
  frailty_variance = 0.2,
  tau = 4
)
aucs <- MCC::CompareAUCs(
  data,
  strata = data$strata,
  tau = 4,
  boot = TRUE,
  perm = TRUE,
  reps = 200,
  alpha = 0.05
)
show(aucs)
```

    ## Marginal Areas:
    ##   arm   n area    se tau
    ## 1   0 100 5.91 0.674   4
    ## 2   1 100 5.06 0.545   4
    ## 
    ## 
    ## CIs:
    ##       method contrast observed    se  lower upper
    ## 1 asymptotic    A1-A0   -0.845 0.867 -2.540 0.854
    ## 3  bootstrap    A1-A0   -0.845 0.865 -2.610 1.000
    ## 2 asymptotic    A1/A0    0.857 0.134  0.630 1.170
    ## 4  bootstrap    A1/A0    0.857 0.137  0.636 1.210
    ## 
    ## 
    ## P-values:
    ##        method contrast observed     p
    ## 1  asymptotic    A1-A0   -0.845 0.330
    ## 3   bootstrap    A1-A0   -0.845 0.348
    ## 5 permutation    A1-A0   -0.845 0.418
    ## 2  asymptotic    A1/A0    0.857 0.325
    ## 4   bootstrap    A1/A0    0.857 0.348
    ## 6 permutation    A1/A0    0.857 0.408

#### Outputs

The `CompareAUCs` object includes:

- `@StratumAreas`: stratum-specific AUCs per arm.

``` r
aucs@StratumAreas
```

    ##   arm strata  n tau     area var_area   se_area strat_weight
    ## 1   0      0 71   4 5.123762 41.65009 0.7659119        0.725
    ## 2   0      1 29   4 7.976865 55.78223 1.3869121        0.275
    ## 3   1      0 74   4 4.860945 29.09823 0.6270720        0.725
    ## 4   1      1 26   4 5.597732 31.10359 1.0937513        0.275

- `@MargAreas`: marginal AUCs per arm.

``` r
aucs@MargAreas
```

    ##   arm   n     area        se tau
    ## 1   0 100 5.908366 0.6736537   4
    ## 2   1 100 5.063561 0.5451197   4

- `@CIs`: confidence intervals for the difference and ratio of AUCs.

``` r
aucs@CIs
```

    ##       method contrast   observed        se      lower     upper
    ## 1 asymptotic    A1-A0 -0.8448043 0.8665822 -2.5432742 0.8536656
    ## 3  bootstrap    A1-A0 -0.8448043 0.8652351 -2.6080192 1.0044493
    ## 2 asymptotic    A1/A0  0.8570156 0.1343891  0.6302478 1.1653760
    ## 4  bootstrap    A1/A0  0.8570156 0.1365854  0.6364443 1.2147199

- `@MCF`: per-arm MCFs.

``` r
head(aucs@MCF)
```

    ##           time         mcf     var_mcf     se_mcf arm
    ## 1 0.0003918607 0.009797297 0.007007054 0.08370814   1
    ## 2 0.0045640871 0.019594595 0.013822133 0.11756757   1
    ## 3 0.0209474391 0.019594595 0.013822133 0.11756757   1
    ## 4 0.0277815206 0.019594595 0.013822133 0.11756757   1
    ## 5 0.0318463015 0.029526101 0.020626585 0.14361958   1
    ## 6 0.0436269205 0.039457608 0.027231065 0.16501838   1

- `@Pvals`: bootstrap and permutation p-values.

``` r
aucs@Pvals
```

    ##        method contrast   observed         p
    ## 1  asymptotic    A1-A0 -0.8448043 0.3296251
    ## 3   bootstrap    A1-A0 -0.8448043 0.3482587
    ## 5 permutation    A1-A0 -0.8448043 0.4179104
    ## 2  asymptotic    A1/A0  0.8570156 0.3251229
    ## 4   bootstrap    A1/A0  0.8570156 0.3482587
    ## 6 permutation    A1/A0  0.8570156 0.4079602

- `@Reps`: bootstrap and permutation test statistics.

#### Influence function

`InfluenceFunction` calculates per-subject influence contributions for
the MCF (`type = "MCF"`) or AUMCF (`type = "AUC"`) at time $\tau$; see
Methods for the definition of $\psi_i(\tau)$.

``` r
psi_auc <- MCC::InfluenceFunction(
  data %>% dplyr::filter(arm == 0),
  tau = 4,
  type = "AUC"
)
head(psi_auc)
```

    ##   idx        psi
    ## 1 101  0.4459027
    ## 2 102 -2.7169950
    ## 3 103 -1.8083700
    ## 4 104 -6.1644698
    ## 5 105 -2.3703085
    ## 6 106  0.3932801

#### Adjusted AUCs

With continuous covariates, `CompareAUCs` uses an augmentation
estimator. Do not supply both `strata` and `covar`; use `model.matrix`
to combine them if needed.

``` r
set.seed(100)
n <- 1000
covariates <- data.frame(
  arm = c(rep(1, n), rep(0, n)),
  x1 = c(stats::rnorm(n, mean = -1), stats::rnorm(n, mean = 1)),
  x2 = c(stats::rnorm(n, mean = 1), stats::rnorm(n, mean = -1))
)
data <- MCC::GenData(
  beta_event = c(log(0.5), log(0.8), log(1.2)),
  covariates = covariates,
  base_death_rate = 0.25,
  base_event_rate = 1,
  frailty_variance = 0.2,
  tau = 4
)
unadj_aucs <- MCC::CompareAUCs(data, tau = 4, alpha = 0.05)
adj_aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  alpha = 0.05,
  covar = data %>% dplyr::select(x1, x2)
)
show(unadj_aucs)
show(adj_aucs)
```

    ## Marginal Areas:
    ##   arm    n area    se tau
    ## 1   0 1000 3.74 0.147   4
    ## 2   1 1000 4.11 0.153   4
    ## 
    ## 
    ## CIs:
    ##       method contrast observed     se   lower upper
    ## 1 asymptotic    A1-A0    0.377 0.2120 -0.0388 0.794
    ## 2 asymptotic    A1/A0    1.100 0.0597  0.9900 1.220
    ## 
    ## 
    ## P-values:
    ##       method contrast observed      p
    ## 1 asymptotic    A1-A0    0.377 0.0755
    ## 2 asymptotic    A1/A0    1.100 0.0757
    ## 
    ## 
    ## Marginal Areas:
    ##   arm    n tau area    se
    ## 1   0 1000   4 3.74 0.147
    ## 2   1 1000   4 4.11 0.153
    ## 
    ## 
    ## CIs:
    ##       method contrast observed    se lower upper
    ## 1 asymptotic    A1-A0    -2.45 0.207 -2.85 -2.04
    ## 
    ## 
    ## P-values:
    ##       method contrast observed        p
    ## 1 asymptotic    A1-A0    -2.45 3.48e-32

#### Plotting

``` r
q_uni <- MCC::PlotMCFs(
  data = data,
  tau = 4,
  x_breaks = seq(0, 4, by = 1),
  color_labs = c("Control", "Treatment")
)
q_uni
```

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

See the [plotting
vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/plotting.pdf)
for additional options.

## Multivariate

Multivariate analysis considers several potentially dependent multiple-
or recurrent-events processes, each subject to right-censoring and
terminal events.

### Methods

For event types $k = 1,\ldots,K$, the arm–process estimand is the AUMCF
$\alpha_{ak}(\tau)$. The multivariate contrast is

$$
\hat\Delta(\tau) = \bigl(\hat\Delta_1(\tau),\ldots,\hat\Delta_K(\tau)\bigr)^\top,
\qquad \hat\Delta_k(\tau) = \hat\alpha_{1k}(\tau) - \hat\alpha_{0k}(\tau).
$$

Optional `process_weights` $w \in \mathbb{R}^K$ defines a secondary
scalar estimand $w^\top\hat\Delta(\tau)$, returned in `@Weighted`.
Per-event jump weighting is not supported for multivariate data; each
recurrent event contributes a jump of size 1.

Arm-specific influence vectors $\hat\phi_{ia} \in \mathbb{R}^K$ are
stacked from per-process influence functions. The arm covariance is

$$
\hat\Sigma_a = \frac{1}{n_a} \sum_{i: A_i=a} \hat\phi_{ia}\hat\phi_{ia}^\top,
$$

and
$\widehat{\mathrm{Cov}}(\hat\Delta) = \hat\Sigma_1/n_1 + \hat\Sigma_0/n_0$.

When patient $i$ is eligible for only a subset of processes, let
$Q_{ik} = I(\text{patient } i \text{ eligible for process } k)$ and
$\hat\pi_{ak} = n_{ak}/n_a$. Per-process influence functions on the
eligible subset are embedded in the arm-level matrix with eligibility
weight $Q_{ik}/\hat\pi_{ak}$; ineligible patients contribute zero.

### Data

Multivariate data are in **long format** with an `event_type` column on
every row, plus `idx`, `arm`, `time`, and `status`. Set
`n_event_types > 1` in `GenData` to simulate such data.

#### Risk sets and terminal events

A subject is at risk for process $k$ if they have at least one row with
non-missing `event_type = k` (recurrent or typed terminal). Subjects
with no type-$k$ rows are excluded from $\hat\Delta_k(\tau)$. Terminal
rows are interpreted as follows:

| `status` | `event_type` | Effect                                           |
|----------|--------------|--------------------------------------------------|
| 0 or 2   | missing      | Ends follow-up for all processes that rely on it |
| 0 or 2   | $k$          | Ends follow-up for process $k$ only              |

A global terminal row does not place a subject in the risk set for
processes without typed rows. To include a subject at risk for $k$ with
zero recurrent type-$k$ events, add a typed censoring or death row with
`event_type = k`.

For each event type $k$ for which a subject is at risk, follow-up must
end with a single observation-terminating event. By default
(`cens_after_last = TRUE`), a censoring event is added after the last
recurrent event for each remaining at-risk event type when no typed
terminal row is present. Set `cens_after_last = FALSE` to leave such
patients at risk indefinitely (with a warning).

### Analyses

#### Multivariate AUCs

``` r
covariates <- data.frame(
  arm = c(rep(1, 100), rep(0, 100))
)
mv_data <- MCC::GenData(
  beta_event = c(log(0.8)),
  covariates = covariates,
  n_event_types = 2L,
  base_event_rate = c(1.0, 0.8),
  death_frailty_variance = 0.2,
  event_frailty_variance = 0.2,
  tau = 4
)
head(mv_data)
```

    ##   idx status      time event_type event_rate arm cens_rate death_rate
    ## 1   1      0 0.9317431         NA         NA   1      0.25  0.1491921
    ## 2   2      1 1.8535908          1  0.1913279   1      0.25  0.1013353
    ## 3   2      2 2.1421942         NA         NA   1      0.25  0.1013353
    ## 4   3      1 1.8915038          1  0.6923577   1      0.25  0.1749469
    ## 5   3      1 2.2586972          1  0.6923577   1      0.25  0.1749469
    ## 6   3      1 2.8287012          2  0.5538862   1      0.25  0.1749469
    ##   death_frailty event_frailty   frailty
    ## 1     0.5967683     1.3802323 0.8236789
    ## 2     0.4053410     0.5900213 0.2391598
    ## 3     0.4053410     0.5900213 0.2391598
    ## 4     0.6997878     1.2367280 0.8654472
    ## 5     0.6997878     1.2367280 0.8654472
    ## 6     0.6997878     1.2367280 0.8654472

``` r
mv_aucs <- MCC::CompareMvAUCs(
  mv_data,
  tau = 4,
  reps = 200
)
show(mv_aucs)
```

    ## Marginal Areas:
    ##    arm event_type  n tau area    se
    ## 1    0          1 62   4 7.96 0.842
    ## 2    0          2 62   4 7.22 0.569
    ## 11   1          1 55   4 6.44 0.547
    ## 21   1          2 55   4 5.94 0.545
    ## 
    ## 
    ## CIs:
    ##        method contrast event_type observed    se lower  upper
    ## 1  asymptotic    A1-A0          1    -1.52 1.000 -3.49 0.4500
    ## 11  bootstrap    A1-A0          1    -1.52 1.110 -3.86 0.5020
    ## 2  asymptotic    A1-A0          2    -1.29 0.788 -2.83 0.2590
    ## 21  bootstrap    A1-A0          2    -1.29 0.775 -2.71 0.0605
    ## 
    ## 
    ## P-values:
    ##        method contrast event_type observed      p
    ## 1  asymptotic    A1-A0          1    -1.52 0.1310
    ## 2  asymptotic    A1-A0          2    -1.29 0.1030
    ## 11  bootstrap    A1-A0          1    -1.52 0.2190
    ## 21  bootstrap    A1-A0          2    -1.29 0.0896

By default, inference is asymptotic (`reps = NULL`); setting `reps` adds
subject-level bootstrap CIs and p-values per event type. Optional
`covars` invoke multivariate augmentation. A weighted scalar contrast is
available via `process_weights`:

``` r
mv_weighted <- MCC::CompareMvAUCs(
  mv_data,
  tau = 4,
  process_weights = c(1, 1)
)
mv_weighted@Weighted
```

    ##   contrast  observed      se     lower    upper          p
    ## 1  w'Delta -2.803201 1.50447 -5.751908 0.145507 0.06242744

#### Plotting

``` r
q_mv <- MCC::PlotMvMCFs(
  mv_data,
  tau = 4,
  x_breaks = seq(0, 4, by = 1),
  event_type_labs = c("Hospitalization", "Heart failure")
)
q_mv
```

![](README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

See the [plotting
vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/plotting.pdf)
for additional multivariate plotting options.
