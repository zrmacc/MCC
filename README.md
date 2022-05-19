# Comparison of mean cumulative count curves via the area under the curve (AUC)

Zachary R. McCaw <br>
Updated: 2022-05-19



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
##   idx status      time arm cens_rate death_rate event_rate   frailty
## 1   1      0 1.7302440   1      0.25  0.2777160  0.8886913 1.1108641
## 2   2      0 1.0019519   1      0.25  0.1287498  0.4119993 0.5149992
## 3   3      2 2.1274945   1      0.25  0.2983774  0.9548077 1.1935096
## 4   4      0 0.3078111   1      0.25  0.5038046  1.6121746 2.0152183
## 5   5      1 0.3701476   1      0.25  0.4659766  1.4911250 1.8639063
## 6   5      1 0.3738888   1      0.25  0.4659766  1.4911250 1.8639063
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

### AUCs

To compare the areas under the mean cumulative count curves up to time $\tau = 4$: 

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
## 1   0 100 5.61 0.561   4
## 2   1 100 5.15 0.622   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -0.460 0.837 -2.100  1.18
## 3  bootstrap    A1-A0   -0.460 0.833 -1.900  1.35
## 2 asymptotic    A1/A0    0.918 0.144  0.675  1.25
## 4  bootstrap    A1/A0    0.918 0.144  0.691  1.25
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.460 0.582
## 3   bootstrap    A1-A0   -0.460 0.657
## 5 permutation    A1-A0   -0.460 0.587
## 2  asymptotic    A1/A0    0.918 0.585
## 4   bootstrap    A1/A0    0.918 0.657
## 6 permutation    A1/A0    0.918 0.587
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
## 1   0 100  5.9 0.713   4
## 2   1 100  4.6 0.514   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0   -1.300 0.879 -3.030 0.418
## 3  bootstrap    A1-A0   -1.300 0.889 -3.050 0.375
## 2 asymptotic    A1/A0    0.779 0.128  0.564 1.080
## 4  bootstrap    A1/A0    0.779 0.130  0.561 1.090
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -1.300 0.138
## 3   bootstrap    A1-A0   -1.300 0.139
## 5 permutation    A1-A0   -1.300 0.129
## 2  asymptotic    A1/A0    0.779 0.129
## 4   bootstrap    A1/A0    0.779 0.139
## 6 permutation    A1/A0    0.779 0.129
```

#### Outputs

The output of `CompareAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

```r
aucs@StratumAreas
```

```
##   arm strata  n tau     area var_area   se_area weight
## 1   0      0 78   4 5.184795 53.84282 0.8308385  0.775
## 2   0      1 22   4 8.372311 40.85371 1.3627130  0.225
## 3   1      0 77   4 4.400773 24.75446 0.5669978  0.775
## 4   1      1 23   4 5.275285 32.25659 1.1842552  0.225
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 5.901986 0.7131739   4
## 2   1 100 4.597538 0.5138992   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##       method contrast   observed        se      lower     upper
## 1 asymptotic    A1-A0 -1.3044477 0.8790389 -3.0273322 0.4184368
## 3  bootstrap    A1-A0 -1.3044477 0.8893716 -3.0468535 0.3749549
## 2 asymptotic    A1/A0  0.7789816 0.1282259  0.5641744 1.0755758
## 4  bootstrap    A1/A0  0.7789816 0.1296619  0.5613733 1.0903481
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##         time        mcf     var_mcf     se_mcf arm
## 1 0.01729211 0.01006494 0.007699022 0.08774407   1
## 2 0.01901668 0.01006494 0.007699022 0.08774407   1
## 3 0.02129013 0.02026230 0.015392641 0.12406708   1
## 4 0.03039660 0.03045967 0.022875552 0.15124666   1
## 5 0.03124162 0.04065704 0.030147753 0.17363108   1
## 6 0.03282630 0.05085441 0.037209246 0.19289698   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##        method contrast   observed         p
## 1  asymptotic    A1-A0 -1.3044477 0.1378228
## 3   bootstrap    A1-A0 -1.3044477 0.1393035
## 5 permutation    A1-A0 -1.3044477 0.1293532
## 2  asymptotic    A1/A0  0.7789816 0.1291764
## 4   bootstrap    A1/A0  0.7789816 0.1393035
## 6 permutation    A1/A0  0.7789816 0.1293532
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
## 1   0 100 5.57 0.561   4
## 2   1 100 5.61 0.626   4
## 
## 
## CIs:
##       method contrast observed    se  lower upper
## 1 asymptotic    A1-A0    0.039 0.840 -1.610  1.69
## 3  bootstrap    A1-A0    0.039 0.795 -1.380  1.57
## 2 asymptotic    A1/A0    1.010 0.151  0.750  1.35
## 4  bootstrap    A1/A0    1.010 0.145  0.776  1.30
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0    0.039 0.963
## 3   bootstrap    A1-A0    0.039 1.000
## 5 permutation    A1-A0    0.039 0.965
## 2  asymptotic    A1/A0    1.010 0.963
## 4   bootstrap    A1/A0    1.010 1.000
## 6 permutation    A1/A0    1.010 0.965
```

### Plotting

The function `PlotMCFs` plots the mean cumulative count curves, comparing two treatment arms. Note that `data` must contain the column `arm`.


```r
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

<img src="README_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />
