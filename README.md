# Comparison of mean cumulative count curves via the area under the curve (AUC)

[![R-CMD-check](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml)

Zachary R. McCaw <br>
Updated: 2026-02-21



### Description

This package provides functions for inference on the difference and ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC curves are estimated using an approach based on the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which accounts for the presence of terminal competing risks. Also see:

* [CICs](https://github.com/zrmacc/CICs) for comparing cumulative incidence curves. 

## Installation


``` r
remotes::install_github("zrmacc/MCC", build_vignettes = TRUE)
```

## Methods

### Estimation of the mean cumulative count curve

For each study arm, the MCC is estimated as follows. Define $N(\cdot)$ as the counting process for events of interest, both terminal and non-terminal, $Y(\cdot)$ as the number of subjects who remain at risk, and $S(\cdot)$ as the probability of not having experienced a terminal event. The MCC $\mu(t)$ at time $t$ is estimated by:

$$
\mu(t) = \int_{0}^{t}\ \hat{S}(u)\ \frac{ dN(u) }{ Y(u) }
$$

Here $\hat{S}(u)$ is the Kaplan-Meier estimate of the probability of being terminal event-free, estimated from *all* terminal events, both those of interest and those regarded as a competing risk; $dN(u)$ is the number of events of interest, both non-terminal and terminal, occurring at time $u$; and $Y(u)$ is the number of subjects who remain at risk, which are subjects who have neither been censored nor experienced a terminal event.

### Standard error calibration

See the [calibration vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/calibration.pdf).

### Data

The function `GenData` simulates example data in the format expected by this package. The recurrent event times are generated from a Poisson process that continues until censoring or death, whichever occurs first. Optionally, a gamma `frailty_variance` may be specified such that the patient-specific event and death rates are correlated. The example data includes 100 patients in each of the treatment and control arms. The maximum duration of follow-up is `tau = 4` (e.g. years). The rate of recurrent events for patients in the treatment arm is 80% the rate for patients in the control arm. 


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

```
##   idx status     time arm cens_rate death_rate event_rate   frailty
## 1   1      2 1.702496   1      0.25 0.33900579  1.0848185 1.3560232
## 2   2      1 2.298376   1      0.25 0.12134480  0.3883034 0.4853792
## 3   2      0 3.251790   1      0.25 0.12134480  0.3883034 0.4853792
## 4   3      1 1.135041   1      0.25 0.06972523  0.2231207 0.2789009
## 5   3      0 4.000000   1      0.25 0.06972523  0.2231207 0.2789009
## 6   4      1 2.415929   1      0.25 0.11716896  0.3749407 0.4686758
```

The essential data are:

* `idx`, the subject index. 
* `time`, the observation time. 
* `status`, coded 0 for censoring, 1 for an event, 2 for death (or any competing terminal event).
* `arm`, coded as 1 for treatment, 0 for reference. 

For analyzing other data sets, arm and status should have the same coding. Each subject should experience an observation-terminating event, i.e. either death or censoring. 

The example data also include:

* `true_death_rate`, the patient-specific terminal event rate, calculated as `frailty` x `base_death_rate` x `exp(covariates %*% beta_death)`. If omitted, `beta_death` is set to zero.
* `true_event_rate`, the patient-specific recurrent event rate, calculated as `frailty` x `base_event_rate` x `exp(covariates %*% beta_event)`. If omitted, `beta_event` is set to zero.
* `frailty`,the patient-specific frailty drawn from a gamma distribution with mean 1 and the specified variance. 


### Observation-terminating events 

In contrast to the time to first event setting, in the multiple or recurrent events setting, a subject may remain at risk after experiencing an event of interest. An *observation-terminating* event, either censoring or the occurrence of a competing risk, is therefore necessary to remove a subject from the risk set. Conversely, a subject who lacks an observation-terminating event is implicitly assumed to remain at risk indefinitely. If a subject *lacks* an observation-terminating event, then by default `CompareAUCs` will add a censoring time immediately after their last event of interest. For example, if the data for subject 1 were:

```
##   idx time status
## 1   1    2      1
## 2   1    3      1
## 3   1    5      1
```
then, for analysis, the subject is assumed to have been censored after the last event, as in the following:

```
##   idx time status
## 1   1    2      1
## 2   1    3      1
## 3   1    5      1
## 4   1    5      0
```

If a subject who lacks an observation-terminating event should, in fact, remain at risk indefinitely, set `cens_after_last = FALSE`.


### Terminal events of interest

Suppose the endpoint of interest includes a fatal event. One such endpoint is heart failure hospitalization (HFH) or cardiovascular (CV)-death. In this setting, it becomes necessary to distinguish non-fatal events of interest (e.g. HFH), after which the subject remains in the risk set, from fatal events of interest (e.g. CV-death), after which the subject is removed from the risk set. To achieve this, a fatal event of interest should be recorded using two records. The first, with `status = 1`, indicates that an event of interest has occurred. The second, with `status = 2`, indicates that the event was terminal. For example, the following data indicate that subject 1 had 3 events of interest, and that the 3rd event, occurring at `time = 5`, was terminal. 

```
##   idx time status
## 1   1    2      1
## 2   1    3      1
## 3   1    5      1
## 4   1    5      2
```

By contrast, the following data indicate that subject 2 had 3 events of interest, none of which was terminal:

```
##   idx time status
## 1   2    2      1
## 2   2    3      1
## 3   2    5      1
```

Note that, by default, subject 2 is assumed to have been censored after their 3rd event of interest, as in the following:

```
##   idx time status
## 1   2    2      1
## 2   2    3      1
## 3   2    5      1
## 4   2    5      0
```

Although censoring (`status = 0`) and a terminal event (`status = 2`) both remove a subject from the risk set, there is an important distinction. Censoring leaves open the possibility that the subject experienced more events of interest in the future, whereas a terminal event precludes the possibility of any future events of interest.


## Analyses

### Single-arm AUC

To calculate the areas under the mean cumulative count curve for a single arm up to time $\tau = 4$:

``` r
auc <- MCC::SingleArmAUC(
  data %>% dplyr::filter(arm == 0),
  boot = TRUE,
  reps = 200,
  tau = 4
)
show(auc)
```

```
## Marginal Areas:
##   arm   n area    se tau
## 1   0 100 5.05 0.621   4
## 
## 
## CIs:
##       method contrast observed    se lower upper
## 1 asymptotic       A0     5.05 0.621  3.83  6.27
## 2  bootstrap       A0     5.05 0.637  3.84  6.35
## 
## 
## P-values:
##       method contrast observed        p
## 1 asymptotic       A0     5.05 4.26e-16
```

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

```
## Marginal Areas:
##   arm   n area    se tau
## 1   0 100 5.05 0.621   4
## 2   1 100 4.21 0.462   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -0.839 0.775 -2.360 0.679
## 3  bootstrap    A1-A0   -0.839 0.771 -2.320 0.504
## 2 asymptotic    A1/A0    0.834 0.137  0.604 1.150
## 4  bootstrap    A1/A0    0.834 0.139  0.615 1.120
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.839 0.279
## 3   bootstrap    A1-A0   -0.839 0.318
## 5 permutation    A1-A0   -0.839 0.308
## 2  asymptotic    A1/A0    0.834 0.271
## 4   bootstrap    A1/A0    0.834 0.318
## 6 permutation    A1/A0    0.834 0.308
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* `boot` indicates to construct bootstrap confidence intervals. 
* `perm` indicates to perform permutation tests for the difference and ratio of AUCs.
* `reps` is the number of simulation replicates. 
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable.
* `alpha` is 1 minus the desired coverage for confidence intervals. 

### Weighted Analysis

Weights may be supplied to control the size of the jump in the cumulative count curve at each event time (i.e. each time with `status == 1`). The following example weights each event by how many events a patient has experienced. For example, if a patient has 3 events before censoring, the first contributes a jump of size 1, the second a jump of size 2, and the third a jump of size 3. Other weighting schemes are of course possible. Note that the weights assigned to censoring (`status == 0`) and terminal event (`status == 2`) records are not used, and may be set to any value. 


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

```
## Visualization of weights for the first 10 records.
## # A tibble: 10 × 4
##      idx  time status weights
##    <dbl> <dbl>  <dbl>   <int>
##  1     1 1.70       2       1
##  2     2 2.30       1       1
##  3     2 3.25       0       2
##  4     3 1.14       1       1
##  5     3 4          0       2
##  6     4 2.42       1       1
##  7     4 3.41       1       2
##  8     4 3.81       0       3
##  9     5 2.23       0       1
## 10     6 0.173      1       1
```


``` r
aucs <- MCC::CompareAUCs(
  data,
  tau = 4,
  alpha = 0.05,
  weights = data$weights
)
show(aucs)
```

```
## Marginal Areas:
##   arm   n  area   se tau
## 1   0 100 12.40 2.56   4
## 2   1 100  8.24 1.44   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -4.170 2.940 -9.930  1.58
## 2 asymptotic    A1/A0    0.664 0.179  0.391  1.13
## 
## 
## P-values:
##       method contrast observed     p
## 1 asymptotic    A1-A0   -4.170 0.155
## 2 asymptotic    A1/A0    0.664 0.129
```

#### Stratified Analysis

`CompareAUCs` also allows for stratified analysis. Consider a data set, similar to that described previously, but with the additional of a binary stratification factor. The event rate for individuals in stratum 1 is increased by 20%.


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

```
## Marginal Areas:
##   arm   n area    se tau
## 1   0 100 5.31 0.572   4
## 2   1 100 4.94 0.510   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -0.371 0.766 -1.870  1.13
## 3  bootstrap    A1-A0   -0.371 0.807 -1.960  1.30
## 2 asymptotic    A1/A0    0.930 0.139  0.694  1.25
## 4  bootstrap    A1/A0    0.930 0.150  0.688  1.29
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.371 0.628
## 3   bootstrap    A1-A0   -0.371 0.627
## 5 permutation    A1-A0   -0.371 0.657
## 2  asymptotic    A1/A0    0.930 0.627
## 4   bootstrap    A1/A0    0.930 0.627
## 6 permutation    A1/A0    0.930 0.657
```

#### Outputs

The output of `CompareAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

``` r
aucs@StratumAreas
```

```
##   arm strata  n tau     area var_area   se_area strat_weight
## 1   0      0 77   4 5.196358 27.22296 0.5945965        0.775
## 2   0      1 23   4 5.684820 52.25091 1.5072420        0.225
## 3   1      0 78   4 4.625223 22.23865 0.5339577        0.775
## 4   1      1 22   4 6.002449 38.54160 1.3235902        0.225
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


``` r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 5.306262 0.5721510   4
## 2   1 100 4.935099 0.5098374   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


``` r
aucs@CIs
```

```
##       method contrast   observed        se      lower    upper
## 1 asymptotic    A1-A0 -0.3711627 0.7663491 -1.8731793 1.130854
## 3  bootstrap    A1-A0 -0.3711627 0.8073335 -1.9563046 1.299140
## 2 asymptotic    A1/A0  0.9300519 0.1388833  0.6940625 1.246281
## 4  bootstrap    A1/A0  0.9300519 0.1502174  0.6876104 1.288097
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


``` r
head(aucs@MCF)
```

```
##          time         mcf     var_mcf     se_mcf arm
## 1 0.004224988 0.009935897 0.007601598 0.08718715   1
## 2 0.010225384 0.009935897 0.007601598 0.08718715   1
## 3 0.010676543 0.009935897 0.007601598 0.08718715   1
## 4 0.017781209 0.019871795 0.015005687 0.12249770   1
## 5 0.018633194 0.019871795 0.015005687 0.12249770   1
## 6 0.030559808 0.029807692 0.022212267 0.14903780   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


``` r
aucs@Pvals
```

```
##        method contrast   observed         p
## 1  asymptotic    A1-A0 -0.3711627 0.6281546
## 3   bootstrap    A1-A0 -0.3711627 0.6268657
## 5 permutation    A1-A0 -0.3711627 0.6567164
## 2  asymptotic    A1/A0  0.9300519 0.6272464
## 4   bootstrap    A1/A0  0.9300519 0.6268657
## 6 permutation    A1/A0  0.9300519 0.6567164
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate continuous covariates. If covariates are provided, then `CompareAUCs` uses an augmentation estimator to adjust for differences between the treatment groups. Note that strata and covariates should not both be provided. If adjustment for both is needed, use `model.matrix` to generate a design matrix including both covariates and stratum indicators, e.g. `model.matrix(~ 0 + covar + strata, data = data)`, then supply the design matrix `covar` argument.


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
##       method contrast observed   se lower upper
## 1 asymptotic    A1-A0   -0.622 0.21 -1.03 -0.21
## 
## 
## P-values:
##       method contrast observed       p
## 1 asymptotic    A1-A0   -0.622 0.00311
```

### Plotting

See the [plotting vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/plotting.pdf).

