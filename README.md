
# Comparison of mean cumulative count curves via the area under the curve (AUC)

[![R-CMD-check](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml)

Zachary R. McCaw <br> Updated: 2026-03-08

### Description

This package provides functions for inference on the difference and
ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC
curves are estimated using an approach based on the method of [Ghosh and
Lin
(2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x),
which accounts for the presence of terminal competing risks. Also see:

- [CICs](https://github.com/zrmacc/CICs) for comparing cumulative
  incidence curves.

## Installation

``` r
remotes::install_github("zrmacc/MCC", build_vignettes = TRUE)
```

## Methods

### Estimation of the mean cumulative count curve

For each study arm, the MCC is estimated as follows. Define $N(\cdot)$
as the counting process for events of interest, both terminal and
non-terminal, $Y(\cdot)$ as the number of subjects who remain at risk,
and $S(\cdot)$ as the probability of not having experienced a terminal
event. The MCC $\mu(t)$ at time $t$ is estimated by:

$$
\mu(t) = \int_{0}^{t}\ \hat{S}(u)\ \frac{ dN(u) }{ Y(u) }
$$

Here $\hat{S}(u)$ is the Kaplan-Meier estimate of the probability of
being terminal event-free, estimated from *all* terminal events, both
those of interest and those regarded as a competing risk; $dN(u)$ is the
number of events of interest, both non-terminal and terminal, occurring
at time $u$; and $Y(u)$ is the number of subjects who remain at risk,
which are subjects who have neither been censored nor experienced a
terminal event.

### Standard error calibration

See the [calibration
vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/calibration.pdf).

### Data

The function `GenData` simulates example data in the format expected by
this package. The recurrent event times are generated from a Poisson
process that continues until censoring or death, whichever occurs first.
Optionally, a gamma `frailty_variance` may be specified such that the
patient-specific event and death rates are correlated. The example data
includes 100 patients in each of the treatment and control arms. The
maximum duration of follow-up is `tau = 4` (e.g. years). The rate of
recurrent events for patients in the treatment arm is 80% the rate for
patients in the control arm.

``` r
library(MCC)
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
    ## 1   1      1 0.3410881   1      0.25  0.1214263  0.3885640 0.4857050
    ## 2   1      0 1.6749694   1      0.25  0.1214263  0.3885640 0.4857050
    ## 3   2      1 2.5421423   1      0.25  0.2278416  0.7290933 0.9113666
    ## 4   2      1 3.1638195   1      0.25  0.2278416  0.7290933 0.9113666
    ## 5   2      0 4.0000000   1      0.25  0.2278416  0.7290933 0.9113666
    ## 6   3      1 0.4335336   1      0.25  0.2988435  0.9562991 1.1953739

The essential data are:

- `idx`, the subject index.
- `time`, the observation time.
- `status`, coded 0 for censoring, 1 for an event, 2 for death (or any
  competing terminal event).
- `arm`, coded as 1 for treatment, 0 for reference.

For analyzing other data sets, arm and status should have the same
coding. Each subject should experience an observation-terminating event,
i.e. either death or censoring.

The example data also include:

- `true_death_rate`, the patient-specific terminal event rate,
  calculated as `frailty` x `base_death_rate` x
  `exp(covariates %*% beta_death)`. If omitted, `beta_death` is set to
  zero.
- `true_event_rate`, the patient-specific recurrent event rate,
  calculated as `frailty` x `base_event_rate` x
  `exp(covariates %*% beta_event)`. If omitted, `beta_event` is set to
  zero.
- `frailty`,the patient-specific frailty drawn from a gamma distribution
  with mean 1 and the specified variance.

### Observation-terminating events

In contrast to the time to first event setting, in the multiple or
recurrent events setting, a subject may remain at risk after
experiencing an event of interest. An *observation-terminating* event,
either censoring or the occurrence of a competing risk, is therefore
necessary to remove a subject from the risk set. Conversely, a subject
who lacks an observation-terminating event is implicitly assumed to
remain at risk indefinitely. If a subject *lacks* an
observation-terminating event, then by default `CompareAUCs` will add a
censoring time immediately after their last event of interest. For
example, if the data for subject 1 were:

    ##   idx time status
    ## 1   1    2      1
    ## 2   1    3      1
    ## 3   1    5      1

then, for analysis, the subject is assumed to have been censored after
the last event, as in the following:

    ##   idx time status
    ## 1   1    2      1
    ## 2   1    3      1
    ## 3   1    5      1
    ## 4   1    5      0

If a subject who lacks an observation-terminating event should, in fact,
remain at risk indefinitely, set `cens_after_last = FALSE`.

### Terminal events of interest

Suppose the endpoint of interest includes a fatal event. One such
endpoint is heart failure hospitalization (HFH) or cardiovascular
(CV)-death. In this setting, it becomes necessary to distinguish
non-fatal events of interest (e.g. HFH), after which the subject remains
in the risk set, from fatal events of interest (e.g. CV-death), after
which the subject is removed from the risk set. To achieve this, a fatal
event of interest should be recorded using two records. The first, with
`status = 1`, indicates that an event of interest has occurred. The
second, with `status = 2`, indicates that the event was terminal. For
example, the following data indicate that subject 1 had 3 events of
interest, and that the 3rd event, occurring at `time = 5`, was terminal.

    ##   idx time status
    ## 1   1    2      1
    ## 2   1    3      1
    ## 3   1    5      1
    ## 4   1    5      2

By contrast, the following data indicate that subject 2 had 3 events of
interest, none of which was terminal:

    ##   idx time status
    ## 1   2    2      1
    ## 2   2    3      1
    ## 3   2    5      1

Note that, by default, subject 2 is assumed to have been censored after
their 3rd event of interest, as in the following:

    ##   idx time status
    ## 1   2    2      1
    ## 2   2    3      1
    ## 3   2    5      1
    ## 4   2    5      0

Although censoring (`status = 0`) and a terminal event (`status = 2`)
both remove a subject from the risk set, there is an important
distinction. Censoring leaves open the possibility that the subject
experienced more events of interest in the future, whereas a terminal
event precludes the possibility of any future events of interest.

## Analyses

### Single-arm AUC

To calculate the areas under the mean cumulative count curve for a
single arm up to time $\tau = 4$:

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

### AUCs

To compare the AUCs of two treatment arms up to time $\tau = 4$:

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

Here:

- `tau` is the truncation time, or the time up to which the AUC is
  calculated.
- `boot` indicates to construct bootstrap confidence intervals.
- `perm` indicates to perform permutation tests for the difference and
  ratio of AUCs.
- `reps` is the number of simulation replicates.
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if
    applicable.
- `alpha` is 1 minus the desired coverage for confidence intervals.

### Weighted Analysis

Weights may be supplied to control the size of the jump in the
cumulative count curve at each event time (i.e. each time with
`status == 1`). The following example weights each event by how many
events a patient has experienced. For example, if a patient has 3 events
before censoring, the first contributes a jump of size 1, the second a
jump of size 2, and the third a jump of size 3. Other weighting schemes
are of course possible. Note that the weights assigned to censoring
(`status == 0`) and terminal event (`status == 2`) records are not used,
and may be set to any value.

``` r
data <- data %>%
  dplyr::group_by(idx) %>%
  dplyr::mutate(weights = dplyr::row_number()) %>%
  dplyr::ungroup()

cat("Visualization of weights for the first 10 records.\n")
data %>%
  dplyr::select(idx, time, status, weights) %>%
  dplyr::slice(1:10)
```

    ## Visualization of weights for the first 10 records.
    ## # A tibble: 10 × 4
    ##      idx  time status weights
    ##    <dbl> <dbl>  <dbl>   <int>
    ##  1     1 0.341      1       1
    ##  2     1 1.67       0       2
    ##  3     2 2.54       1       1
    ##  4     2 3.16       1       2
    ##  5     2 4          0       3
    ##  6     3 0.434      1       1
    ##  7     3 0.571      1       2
    ##  8     3 1.49       1       3
    ##  9     3 1.86       0       4
    ## 10     4 0.597      2       1

``` r
aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  alpha = 0.05,
  weights = data$weights
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

#### Stratified Analysis

`CompareAUCs` also allows for stratified analysis. Consider a data set,
similar to that described previously, but with the additional of a
binary stratification factor. The event rate for individuals in stratum
1 is increased by 20%.

``` r
# Generate data with strata.
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

# Stratified AUC analysis.
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

The output of `CompareAUCs` is an object with these slots.

- `@StratumAreas` containing the stratum-specific AUCs for each arm.

``` r
aucs@StratumAreas
```

    ##   arm strata  n tau     area var_area   se_area strat_weight
    ## 1   0      0 71   4 5.123762 41.65009 0.7659119        0.725
    ## 2   0      1 29   4 7.976865 55.78223 1.3869121        0.275
    ## 3   1      0 74   4 4.860945 29.09823 0.6270720        0.725
    ## 4   1      1 26   4 5.597732 31.10359 1.0937513        0.275

- `@MargAreas` containing the AUCs for each arm, marginalized over any
  strata.

``` r
aucs@MargAreas
```

    ##   arm   n     area        se tau
    ## 1   0 100 5.908366 0.6736537   4
    ## 2   1 100 5.063561 0.5451197   4

- `@CIs` containing confidence intervals for the difference and ratio of
  AUCs.

``` r
aucs@CIs
```

    ##       method contrast   observed        se      lower     upper
    ## 1 asymptotic    A1-A0 -0.8448043 0.8665822 -2.5432742 0.8536656
    ## 3  bootstrap    A1-A0 -0.8448043 0.8652351 -2.6080192 1.0044493
    ## 2 asymptotic    A1/A0  0.8570156 0.1343891  0.6302478 1.1653760
    ## 4  bootstrap    A1/A0  0.8570156 0.1365854  0.6364443 1.2147199

- `@MCF` containing the per arm mean cumulative count curve, averaged
  across strata.

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

- `@Pvals` containing the bootstrap and permutation p-values.

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

- `@Reps` is a list containing the bootstrap and permutation test
  statistics.

### Influence function

Let $\theta(\tau)$ denote the MCF or AUMCF at time $\tau$. The influence
function $\psi_{i}(\tau)$ measures the contribution of subject $i$ to
the estimation error:

$$
\hat{\theta}(\tau) - \theta(\tau) = \frac{1}{n}\sum_{i=1}^{n}\psi_{i}(\tau)
$$

The influence function contributions are useful for estimating the
variance of $\hat{\theta}(\tau)$:

$$
\hat{\mathbb{V}}(\hat{\theta}) = \frac{1}{n^2}\sum_{i=1}^{n}\psi_{i}^{2}(\tau).
$$

`InfluenceFunction` can be used to calculate the influence contributions
of each subject to the MCF (`type = "MCF"`) or AUMCF (`type = "AUC"`) at
a given time $\tau$.

``` r
# Influence for AUC up to tau (one row per subject).
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

``` r
# Influence for MCF at tau.
psi_mcf <- MCC::InfluenceFunction(
  data %>% dplyr::filter(arm == 0),
  tau = 4,
  type = "MCF"
)
head(psi_mcf)
```

    ##   idx         psi
    ## 1 101  0.03212038
    ## 2 102 -1.76936792
    ## 3 103  1.48895373
    ## 4 104 -2.58604980
    ## 5 105 -0.39307438
    ## 6 106 -0.86405517

### Adjusted AUCs

The previous estimator allows for stratification, but a different
approach is needed to accommodate continuous covariates. If covariates
are provided, then `CompareAUCs` uses an augmentation estimator to
adjust for differences between the treatment groups. Note that strata
and covariates should not both be provided. If adjustment for both is
needed, use `model.matrix` to generate a design matrix including both
covariates and stratum indicators,
e.g. `model.matrix(~ 0 + covar + strata, data = data)`, then supply the
design matrix `covar` argument.

``` r
set.seed(100)

# Generate data with a continuous covariate.
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

# Unadjusted.
paste("Unadjusted AUCs:")
unadj_aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  alpha = 0.05
)
show(unadj_aucs)

# Adjusted.
paste("Adjusted AUCs:")
adj_aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  alpha = 0.05,
  covar = data %>% dplyr::select(x1, x2)
)
show(adj_aucs)
```

    ## [1] "Unadjusted AUCs:"
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
    ## [1] "Adjusted AUCs:"
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

### Plotting

See the [plotting
vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/plotting.pdf).
