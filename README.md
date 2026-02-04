# Comparison of mean cumulative count curves via the area under the curve (AUC)

Zachary R. McCaw <br>
Updated: 2026-02-04

<!-- badges: start -->
[![R-CMD-check](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/MCC/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->



### Description

This package provides functions for inference on the difference and ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC curves are estimated using an approach based on the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which accounts for the presence of terminal competing risks. Also see:

* [CICs](https://github.com/zrmacc/CICs) for comparing cumulative incidence curves. 

## Installation


``` r
devtools::install_github(repo = 'zrmacc/MCC')
```

## Methods

### Estimation of the mean cumulative count curve

For each study arm, the MCC is estimated as follows. Define $N(\cdot)$ as the counting process for events of interest, both terminal and non-terminal, $Y(\cdot)$ as the number of subjects who remain at risk, and $S(\cdot)$ as the probability of not having experienced a terminal event. The MCC $\mu(t)$ at time $t$ is estimated by:
$$
\mu(t) = \int_{0}^{t} \hat{S}(u) \frac{ dN(u) }{ Y(u) }
$$

Here $\hat{S}(u)$ is the Kaplan-Meier estimate of the probability of being terminal event-free, estimated from *all* terminal events, both those of interest and those regarded as a competing risk; $dN(u)$ is the number of events of interest, both non-terminal and terminal, occurring at time $u$; and $Y(u)$ is the number of subjects who remain at risk, which are subjects who have neither been censored nor experienced a terminal event.




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
##   idx status     time arm cens_rate death_rate event_rate  frailty
## 1   1      1 2.018025   1      0.25  0.3429874  1.0975598 1.371950
## 2   1      1 2.686707   1      0.25  0.3429874  1.0975598 1.371950
## 3   1      1 2.804273   1      0.25  0.3429874  1.0975598 1.371950
## 4   1      0 2.848668   1      0.25  0.3429874  1.0975598 1.371950
## 5   2      1 3.002166   1      0.25  0.2287110  0.7318752 0.914844
## 6   2      2 3.355392   1      0.25  0.2287110  0.7318752 0.914844
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
## 1   0 100 6.06 0.581   4
## 
## 
## CIs:
##       method contrast observed    se lower upper
## 1 asymptotic       A0     6.06 0.581  4.92  7.19
## 2  bootstrap       A0     6.06 0.562  4.97  7.10
## 
## 
## P-values:
##       method contrast observed        p
## 1 asymptotic       A0     6.06 1.94e-25
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
## 1   0 100 6.06 0.581   4
## 2   1 100 4.69 0.443   4
## 
## 
## CIs:
##       method contrast observed     se  lower   upper
## 1 asymptotic    A1-A0   -1.360 0.7310 -2.800  0.0671
## 3  bootstrap    A1-A0   -1.360 0.6630 -2.640 -0.1050
## 2 asymptotic    A1/A0    0.775 0.1040  0.595  1.0100
## 4  bootstrap    A1/A0    0.775 0.0956  0.611  0.9790
## 
## 
## P-values:
##        method contrast observed      p
## 1  asymptotic    A1-A0   -1.360 0.0617
## 3   bootstrap    A1-A0   -1.360 0.0398
## 5 permutation    A1-A0   -1.360 0.0896
## 2  asymptotic    A1/A0    0.775 0.0578
## 4   bootstrap    A1/A0    0.775 0.0398
## 6 permutation    A1/A0    0.775 0.0896
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* `boot` indicates to construct bootstrap confidence intervals. 
* `perm` indicates to perform permutation tests for the difference and ratio of AUCs.
* `reps` is the number of simulation replicates. 
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable.
* `alpha` is 1 minus the desired coverage for confidence intervals. 

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
## 1   0 100 5.89 0.622   4
## 2   1 100 5.07 0.467   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -0.817 0.778 -2.340 0.708
## 3  bootstrap    A1-A0   -0.817 0.750 -2.210 0.717
## 2 asymptotic    A1/A0    0.861 0.121  0.654 1.130
## 4  bootstrap    A1/A0    0.861 0.120  0.667 1.140
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.817 0.294
## 3   bootstrap    A1-A0   -0.817 0.328
## 5 permutation    A1-A0   -0.817 0.308
## 2  asymptotic    A1/A0    0.861 0.287
## 4   bootstrap    A1/A0    0.861 0.328
## 6 permutation    A1/A0    0.861 0.308
```

#### Outputs

The output of `CompareAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

``` r
aucs@StratumAreas
```

```
##   arm strata  n tau     area var_area   se_area strat_weight
## 1   0      0 80   4 5.576907 34.57962 0.6574536         0.79
## 2   0      1 20   4 7.073537 53.29044 1.6323363         0.21
## 3   1      0 78   4 4.843003 18.90140 0.4922658         0.79
## 4   1      1 22   4 5.943331 33.55706 1.2350388         0.21
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


``` r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 5.891199 0.6223099   4
## 2   1 100 5.074072 0.4674420   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


``` r
aucs@CIs
```

```
##       method contrast   observed        se      lower    upper
## 1 asymptotic    A1-A0 -0.8171272 0.7783134 -2.3425934 0.708339
## 3  bootstrap    A1-A0 -0.8171272 0.7503167 -2.2090991 0.717044
## 2 asymptotic    A1/A0  0.8612970 0.1207208  0.6544058 1.133597
## 4  bootstrap    A1/A0  0.8612970 0.1199138  0.6673856 1.144192
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


``` r
head(aucs@MCF)
```

```
##          time        mcf     var_mcf     se_mcf arm
## 1 0.008839499 0.00000000 0.000000000 0.00000000   1
## 2 0.014325674 0.01012821 0.007898702 0.08887464   1
## 3 0.025827966 0.01967366 0.009811746 0.09905426   1
## 4 0.026231461 0.02980186 0.017505286 0.13230754   1
## 5 0.047070046 0.03993007 0.024993666 0.15809385   1
## 6 0.047218654 0.03993007 0.024993666 0.15809385   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


``` r
aucs@Pvals
```

```
##        method contrast   observed         p
## 1  asymptotic    A1-A0 -0.8171272 0.2937783
## 3   bootstrap    A1-A0 -0.8171272 0.3283582
## 5 permutation    A1-A0 -0.8171272 0.3084577
## 2  asymptotic    A1/A0  0.8612970 0.2867345
## 4   bootstrap    A1/A0  0.8612970 0.3283582
## 6 permutation    A1/A0  0.8612970 0.3084577
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate continuous covariates. If covariates are provided, then `CompareAUCs` uses an augmentation estimator to adjust for differences between the treatment groups. Note that strata and covariates should not both be provided. If adjustment for both is needed, use `model.matrix` to generate a design matrix including both covariates and stratum indicators, e.g. `model.matrix(~ 0 + covar + strata, data = data)`, then supply the design matrix `covar` argument.


``` r
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
## 1   0 100 6.13 0.671   4
## 2   1 100 4.81 0.544   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -1.320 0.864 -3.020 0.368
## 3  bootstrap    A1-A0   -1.320 0.874 -3.080 0.573
## 2 asymptotic    A1/A0    0.784 0.123  0.576 1.070
## 4  bootstrap    A1/A0    0.784 0.126  0.577 1.110
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -1.320 0.125
## 3   bootstrap    A1-A0   -1.320 0.129
## 5 permutation    A1-A0   -1.320 0.159
## 2  asymptotic    A1/A0    0.784 0.122
## 4   bootstrap    A1/A0    0.784 0.129
## 6 permutation    A1/A0    0.784 0.159
```

### Plotting

The function `PlotMCFs` plots the mean cumulative count curves, comparing two treatment arms. Note that `data` must contain the column `arm`.


``` r
q <- MCC::PlotMCFs(data)
q_nar <- MCC::PlotNARs(
  data = data,
  x_breaks = seq(from = 0, to = 4)
)
q_main <- cowplot::plot_grid(
  plotlist = list(q, q_nar),
  nrow = 2,
  align = "v",
  axis = "l",
  rel_heights = c(3, 1)
)
show(q_main)
```

<img src="README_files/figure-html/unnamed-chunk-17-1.png" alt="" style="display: block; margin: auto;" />
