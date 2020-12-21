# Compare Mean Cumulative Count Curves

Zachary McCaw <br>
Updated: 2020-12-21



### Description

This package provides functions for inference on the difference and ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC curves are estimated using the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which allows for the occurrence of terminal events such as death. Also see:

* [CICs](https://github.com/zrmacc/CICs) for comparing cumulative incidence curves. 

## Installation


```r
devtools::install_github(repo = 'zrmacc/MCC')
```

## Examples

### Data

The function `GenData` simulates example data in the format expected by this package. The censoring, death, and recurrent event gap times are drawn from independent exponential distributions. The example data includes 100 patients in each of the treatment and control arms. Observation of a given patient continues until censoring, death, or time `tau = 4`, whichever occurs first. The rate of recurrent events for patients in the treatment arm is 80% the rate for patients in the control arm. 


```r
library(MCC)
data <- MCC::GenData(
  n1 = 100,
  n0 = 100,
  treatment_effect = log(0.8),
  tau = 4
)
head(data)
```

```
##   idx       time status arm      covar strata true_event_rate
## 1   1 1.24880530      2   1 -0.1191161      0       0.7790161
## 2   2 0.92041765      1   1  1.7388641      0       1.1792431
## 3   2 1.20365788      0   1  1.7388641      0       1.1792431
## 4   3 0.04608151      1   1 -0.2368181      0       0.7588220
## 5   3 1.27370805      1   1 -0.2368181      0       0.7588220
## 6   3 1.56318526      1   1 -0.2368181      0       0.7588220
```

The essential data are:

* `idx`, the subject index. 
* `time`, the observation time. 
* `status`, coded 0 for censoring, 1 for an event, 2 for death (or any competing terminal event).
* `arm`, coded as 1 for treatment, 0 for reference. 

For analyzing other data sets, arm and status should have the same coding. Each subject should experience an observation-terminating event, i.e. either death or censoring. If the last appearance of a subject in the data set has status 1, then a censoring time is added immediately after this recurrence. For example, if the data for subject 1 is:

```
##   idx time status
## 1   1    2      1
## 2   1    3      1
## 3   1    5      1
```
then, for analysis, the subject's is assumed to be censored after the last event:

```
##   idx time status
## 1   1    2      1
## 2   1    3      1
## 3   1    5      1
## 4   1    5      0
```

If instead the last recurrence is fatal, encode the input data as:

```
##   idx time status
## 1   1    2      1
## 2   1    3      1
## 3   1    5      1
## 4   1    5      2
```

The example data also include:

* `covar`, a standard normal covariate that decreases the event rate. 
* `strat`, a two-level stratification factor that increases the event rate.
* `true_event_rate`, the patient-specific recurrent event rate.

### AUCs

To compare the areas under the mean cumulative count curves up to time $\tau = 4$: 

```r
aucs <- MCC::CompareAUCs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  idx = data$idx,
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
## 1   0 100 6.72 0.582   4
## 2   1 100 4.72 0.462   4
## 
## 
## CIs:
##       method contrast observed     se  lower  upper
## 1 asymptotic    A1-A0   -2.000 0.7430 -3.450 -0.539
## 3  bootstrap    A1-A0   -2.000 0.7950 -3.880 -0.662
## 2 asymptotic    A1/A0    0.703 0.0919  0.544  0.908
## 4  bootstrap    A1/A0    0.703 0.0973  0.519  0.890
## 
## 
## P-values:
##        method contrast observed       p
## 1  asymptotic    A1-A0   -2.000 0.00724
## 3   bootstrap    A1-A0   -2.000 0.00995
## 5 permutation    A1-A0   -2.000 0.01990
## 2  asymptotic    A1/A0    0.703 0.00699
## 4   bootstrap    A1/A0    0.703 0.00995
## 6 permutation    A1/A0    0.703 0.01990
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* `boot` indicates to construct bootstrap confidence intervals. 
* `perm` indicates to perform permutation tests for the difference and ratio of AUCs.
* `reps` is the number of simulation replicates. 
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable.
* `alpha` is 1 minus the desired coverage for confidence intervals. 

Note that the `strata` argument may be omitted for unstratified data. 

#### Outputs

The output of `CompareAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

```r
aucs@StratumAreas
```

```
##   arm strata  n tau     area var_area   se_area weight
## 1   0      0 82   4 6.521333 32.07541 0.6254307  0.815
## 2   0      1 18   4 7.584617 41.42754 1.5170795  0.185
## 3   1      0 81   4 4.806613 22.18537 0.5233483  0.815
## 4   1      1 19   4 4.352797 17.50931 0.9599702  0.185
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 6.718040 0.5818853   4
## 2   1 100 4.722657 0.4620246   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##       method contrast   observed         se      lower      upper
## 1 asymptotic    A1-A0 -1.9953835 0.74300550 -3.4516475 -0.5391195
## 3  bootstrap    A1-A0 -1.9953835 0.79473518 -3.8838068 -0.6617801
## 2 asymptotic    A1/A0  0.7029813 0.09185471  0.5441542  0.9081668
## 4  bootstrap    A1/A0  0.7029813 0.09731066  0.5190986  0.8901389
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##          time       mcf     var_mcf     se_mcf arm
## 1 0.004062820 0.0000000 0.000000000 0.00000000   1
## 2 0.005559498 0.0000000 0.000000000 0.00000000   1
## 3 0.008942332 0.0101875 0.008301515 0.09111265   1
## 4 0.019738298 0.0203750 0.016392865 0.12803463   1
## 5 0.025924152 0.0305625 0.024274051 0.15580132   1
## 6 0.027934783 0.0305625 0.024274051 0.15580132   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##        method contrast   observed           p
## 1  asymptotic    A1-A0 -1.9953835 0.007240906
## 3   bootstrap    A1-A0 -1.9953835 0.009950249
## 5 permutation    A1-A0 -1.9953835 0.019900498
## 2  asymptotic    A1/A0  0.7029813 0.006993068
## 4   bootstrap    A1/A0  0.7029813 0.009950249
## 6 permutation    A1/A0  0.7029813 0.019900498
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate continuous covariates. If covariates are provided, then `CompareAUCs` uses an augmentation estimator to adjust for differences between the treatment groups. Note that strata and covariates should not both be provided. If adjustment for both is needed, use `model.matrix` to generate a design matrix including both covariates and stratum indicators, e.g. `model.matrix(~ 0 + covar + strata, data = data)`, then supply the design matrix `covar` argument.


```r
aucs <- MCC::CompareAUCs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  idx = data$idx,
  covar = data$covar,
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
##   arm   n tau area    se
## 1   0 100   4 6.72 0.586
## 2   1 100   4 4.75 0.467
## 
## 
## CIs:
##       method contrast observed    se lower   upper
## 1 asymptotic    A1-A0    -1.01 0.529 -2.05  0.0235
## 2  bootstrap    A1-A0    -1.01 0.532 -2.32 -0.0756
## 
## 
## P-values:
##        method contrast observed      p
## 1  asymptotic    A1-A0    -1.01 0.0554
## 2   bootstrap    A1-A0    -1.01 0.0597
## 3 permutation    A1-A0    -1.01 0.1090
```
