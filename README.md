# Comparison of mean cumulative count curves via the area under the curve (AUC)

Zachary R. McCaw <br>
Updated: 2024-02-24



### Description

This package provides functions for inference on the difference and ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC curves are estimated using an approach based on the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which accounts for the presence of terminal competing risks. Also see:

* [CICs](https://github.com/zrmacc/CICs) for comparing cumulative incidence curves. 

## Installation


```r
devtools::install_github(repo = 'zrmacc/MCC')
```

## Methods

### Estimation of the mean cumulative count curve

For each study arm, the MCC is estimated as follows. Define $N(\cdot)$ as the counting process for events of interest, both terminal and non-terminal, $Y(\cdot)$ as the number of subjects who remain at risk, and $S(\cdot)$ as the probability of not having experienced a terminal event. The MCC $\mu(t)$ at time $t$ is estimated by:
$$
\mu(t) = \int_{0}^{t} \hat{S}(u) \frac{ dN(u) }{ Y(u) }
$$

Here $\hat{S}(u)$ is the Kaplan-Meier estimate of the probability of being terminal event-free, estimated from *all* terminal events, both those of interest and those regarded as a competing risk; $dN(u)$ is the number of events of interest, both non-terminal and terminal, occurring at time $u$; and $Y(u)$ is the number of subjects who remain at risk, which are subjects who have neither been censored nor experienced a terminal event.

### Standard error calibration

See the [calibration vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/calibration.pdf).

### Data

The function `GenData` simulates example data in the format expected by this package. The recurrent event times are generated from a Poisson process that continues until censoring or death, whichever occurs first. Optionally, a gamma `frailty_variance` may be specified such that the patient-specific event and death rates are correlated. The example data includes 100 patients in each of the treatment and control arms. The maximum duration of follow-up is `tau = 4` (e.g. years). The rate of recurrent events for patients in the treatment arm is 80% the rate for patients in the control arm. 


```r
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
##   idx status       time arm cens_rate death_rate event_rate   frailty
## 1   1      1 0.08845737   1      0.25  0.1888051  0.6041762 0.7552202
## 2   1      2 0.23565323   1      0.25  0.1888051  0.6041762 0.7552202
## 3   2      1 0.16410134   1      0.25  0.1884027  0.6028886 0.7536108
## 4   2      1 0.66910083   1      0.25  0.1884027  0.6028886 0.7536108
## 5   2      2 0.87297407   1      0.25  0.1884027  0.6028886 0.7536108
## 6   3      2 0.02295830   1      0.25  0.3140817  1.0050614 1.2563267
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

```r
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
## 1   0 100 6.29 0.568   4
## 
## 
## CIs:
##       method contrast observed    se lower upper
## 1 asymptotic       A0     6.29 0.568  5.18  7.40
## 2  bootstrap       A0     6.29 0.641  5.31  7.65
## 
## 
## P-values:
##       method contrast observed        p
## 1 asymptotic       A0     6.29 1.55e-28
```

### AUCs

To compare the AUCs of two treatment arms up to time $\tau = 4$: 

```r
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
## 1   0 100 6.29 0.568   4
## 2   1 100 4.31 0.476   4
## 
## 
## CIs:
##       method contrast observed     se  lower  upper
## 1 asymptotic    A1-A0   -1.980 0.7410 -3.430 -0.530
## 3  bootstrap    A1-A0   -1.980 0.7120 -3.350 -0.616
## 2 asymptotic    A1/A0    0.685 0.0978  0.518  0.906
## 4  bootstrap    A1/A0    0.685 0.0980  0.525  0.881
## 
## 
## P-values:
##        method contrast observed       p
## 1  asymptotic    A1-A0   -1.980 0.00748
## 3   bootstrap    A1-A0   -1.980 0.01990
## 5 permutation    A1-A0   -1.980 0.00995
## 2  asymptotic    A1/A0    0.685 0.00801
## 4   bootstrap    A1/A0    0.685 0.01990
## 6 permutation    A1/A0    0.685 0.00995
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* `boot` indicates to construct bootstrap confidence intervals. 
* `perm` indicates to perform permutation tests for the difference and ratio of AUCs.
* `reps` is the number of simulation replicates. 
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable.
* `alpha` is 1 minus the desired coverage for confidence intervals. 

### Weighted Analysis

Weights may be supplied to control the size of the jump in the cumulative count curve at each event time (i.e. each time with `status == 1`). The following example weights each event by how many events a patient has experienced. For example, if a patient has 3 events before censoring, the first contributes a jump of size 1, the second of size 2, and the third of size 3. Other weighting schemes are of course possible. Note that the weights assigned to censoring (`status == 0`) and terminal event (`status == 2`) records are not used, and may be set to any value.


```r
data <- data %>%
  dplyr::group_by(idx) %>%
  dplyr::mutate(weights = dplyr::row_number()) %>%
  dplyr::ungroup()

# Example.
data %>%
  dplyr::select(idx, time, status, weights) %>%
  dplyr::slice(1:10)
```

```
## # A tibble: 10 × 4
##      idx   time status weights
##    <dbl>  <dbl>  <dbl>   <int>
##  1     1 0.0885      1       1
##  2     1 0.236       2       2
##  3     2 0.164       1       1
##  4     2 0.669       1       2
##  5     2 0.873       2       3
##  6     3 0.0230      2       1
##  7     4 0.0975      0       1
##  8     5 0.257       0       1
##  9     6 0.0402      1       1
## 10     6 1.29        1       2
```


```r
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
## 1   0 100 13.70 1.97   4
## 2   1 100  8.42 1.37   4
## 
## 
## CIs:
##       method contrast observed    se  lower  upper
## 1 asymptotic    A1-A0   -5.250 2.400 -9.960 -0.541
## 2 asymptotic    A1/A0    0.616 0.134  0.402  0.944
## 
## 
## P-values:
##       method contrast observed      p
## 1 asymptotic    A1-A0   -5.250 0.0289
## 2 asymptotic    A1/A0    0.616 0.0260
```

#### Stratified Analysis

`CompareAUCs` also allows for stratified analysis. Consider a data set, similar to that described previously, but with the additional of a binary stratification factor. The event rate for individuals in stratum 1 is increased by 20%.


```r
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
## 1   0 100 4.53 0.505   4
## 2   1 100 5.22 0.497   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0    0.694 0.709 -0.695  2.08
## 3  bootstrap    A1-A0    0.694 0.696 -0.877  1.89
## 2 asymptotic    A1/A0    1.150 0.169  0.865  1.54
## 4  bootstrap    A1/A0    1.150 0.163  0.819  1.47
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0    0.694 0.327
## 3   bootstrap    A1-A0    0.694 0.348
## 5 permutation    A1-A0    0.694 0.318
## 2  asymptotic    A1/A0    1.150 0.331
## 4   bootstrap    A1/A0    1.150 0.348
## 6 permutation    A1/A0    1.150 0.318
```

#### Outputs

The output of `CompareAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

```r
aucs@StratumAreas
```

```
##   arm strata  n tau     area var_area   se_area strat_weight
## 1   0      0 79   4 4.475710 26.88668 0.5833847         0.81
## 2   0      1 21   4 4.762774 18.28877 0.9332170         0.19
## 3   1      0 83   4 4.629344 22.07835 0.5157560         0.81
## 4   1      1 17   4 7.760514 34.30218 1.4204843         0.19
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 4.530252 0.5047126   4
## 2   1 100 5.224266 0.4973601   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##       method contrast  observed        se      lower    upper
## 1 asymptotic    A1-A0 0.6940141 0.7085915 -0.6947997 2.082828
## 3  bootstrap    A1-A0 0.6940141 0.6960016 -0.8768036 1.890103
## 2 asymptotic    A1/A0 1.1531955 0.1689951  0.8652937 1.536888
## 4  bootstrap    A1/A0 1.1531955 0.1630228  0.8191117 1.473231
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##           time        mcf     var_mcf     se_mcf arm
## 1 2.029460e-06 0.01117647 0.001998616 0.04470588   1
## 2 3.797853e-04 0.02235294 0.003747405 0.06121605   1
## 3 3.983363e-03 0.03352941 0.005246367 0.07243181   1
## 4 1.167436e-02 0.04328845 0.013055947 0.11426262   1
## 5 1.297160e-02 0.04328845 0.013055947 0.11426262   1
## 6 1.868159e-02 0.05316650 0.020860825 0.14443277   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##        method contrast  observed         p
## 1  asymptotic    A1-A0 0.6940141 0.3273687
## 3   bootstrap    A1-A0 0.6940141 0.3482587
## 5 permutation    A1-A0 0.6940141 0.3184080
## 2  asymptotic    A1/A0 1.1531955 0.3307283
## 4   bootstrap    A1/A0 1.1531955 0.3482587
## 6 permutation    A1/A0 1.1531955 0.3184080
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate continuous covariates. If covariates are provided, then `CompareAUCs` uses an augmentation estimator to adjust for differences between the treatment groups. Note that strata and covariates should not both be provided. If adjustment for both is needed, use `model.matrix` to generate a design matrix including both covariates and stratum indicators, e.g. `model.matrix(~ 0 + covar + strata, data = data)`, then supply the design matrix `covar` argument.


```r
# Generate data with a continuous covariate.
covariates <- data.frame(
  arm = c(rep(1, 100), rep(0, 100)),
  covar = stats::rnorm(200)
)
data <- MCC::GenData(
  beta_event = c(log(0.8), log(1.2)),
  covariates = covariates,
  frailty_variance = 0.2,
  tau = 4
)

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
## 1   0 100 5.32 0.600   4
## 2   1 100 5.11 0.674   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -0.215 0.903 -1.980  1.56
## 3  bootstrap    A1-A0   -0.215 0.806 -1.640  1.58
## 2 asymptotic    A1/A0    0.960 0.167  0.683  1.35
## 4  bootstrap    A1/A0    0.960 0.147  0.710  1.32
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.215 0.812
## 3   bootstrap    A1-A0   -0.215 0.726
## 5 permutation    A1-A0   -0.215 0.726
## 2  asymptotic    A1/A0    0.960 0.812
## 4   bootstrap    A1/A0    0.960 0.726
## 6 permutation    A1/A0    0.960 0.726
```

### Plotting

See the [plotting vignette](https://github.com/zrmacc/MCC/blob/master/vignettes/plotting.pdf).

