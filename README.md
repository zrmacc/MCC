# Compare Mean Cumulative Count Curves

Zachary McCaw <br>
Updated: 2020-12-17



### Description

This package provides functions for inference on the difference and ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC curves are estimated using the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which allows for the occurrence of terminal events such as death. Also see:
* [CICs](https://github.com/zrmacc/CICs) for comparing cumulative incidence curves. 

## Installation


```r
devtools::install_github(repo = 'zrmacc/MCC')
```

## Examples

### Data

The function `GenData` simulates example data in the format expected by this package. The censoring, death, and recurrent event gap times are drawn from independent exponential distributions. The example data includes 100 patients in each of the treatment and control arms. Observation of a given patient continues until censoring, death, or time `tau = 10`, whichever occurs first. The rate of recurrent events for patients in the treatment arm is 80% the rate for patients in the control arm. 


```r
library(MCC)
data <- MCC::GenData(
  n1 = 100,
  n0 = 100,
  treatment_effect = log(0.8),
  tau = 10
)
head(data)
```

```
##   idx      time status arm     covar strata true_event_rate
## 1   1 1.3310668      0   1 -2.148078      0       0.4953587
## 2   2 0.2454117      1   1  0.724350      0       0.9403440
## 3   2 0.7317782      2   1  0.724350      0       0.9403440
## 4   3 0.9460894      2   1  1.131180      0       1.0297045
## 5   4 1.2020167      1   1 -1.509513      0       0.5712195
## 6   4 3.7312608      1   1 -1.509513      0       0.5712195
```

The essential data are:

* `idx`, the subject index. 
* `time`, the observation time. 
* `status`, coded 0 for censoring, 1 for an event, 2 for death.
* `arm`, coded as 1 for treatment, 0 for reference. 

For analyses of other data sets, arm and status should have the same coding. Each subject should experience at most one of censoring or death. Subjects who experience neither censoring nor death are assumed to remain at risk throughout follow-up. 

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
  reps = 100,
  alpha = 0.05
)
show(aucs)
```

```
## Marginal Areas:
##   arm   n area    se tau
## 1   0 100 6.44 0.523   4
## 2   1 100 5.07 0.563   4
## 
## 
## CIs:
##          method contrast observed     se  lower   upper
## 1    asymptotic    A1-A0   -1.360 0.7690 -2.870  0.1420
## 11    bootstrap    A1-A0   -1.360 0.4040 -1.500 -0.0501
## 2    asymptotic    A1/A0    0.788 0.1080  0.602  1.0300
## 2.5%  bootstrap    A1/A0    0.788 0.0979  0.632  0.9840
## 
## 
## P-values:
##        method contrast observed      p
## 1  asymptotic    A1-A0   -1.360 0.0758
## 3   bootstrap    A1-A0   -1.360 0.0396
## 5 permutation    A1-A0   -1.360 0.0792
## 2  asymptotic    A1/A0    0.788 0.0833
## 4   bootstrap    A1/A0    0.788 0.0396
## 6 permutation    A1/A0    0.788 0.0792
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* `boot` indicates to construct bootstrap confidence intervals. 
* `perm` indicates to perform permutation tests for the difference and ratio of AUCs.
* `reps` is the number of simulation replicates. 
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable.
* `alpha` is 1 minus the desired coverage for confidence intervals. 
* `seed` allows for reproducible bootstrap/permutation replicates.

Note that the `strata` argument may be omitted for unstratified data. 

#### Outputs

The output of `CompareAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

```r
aucs@StratumAreas
```

```
##   arm strata  n tau     area var_area   se_area weight
## 1   0      0 76   4 5.775588 24.20571 0.5643546   0.77
## 2   0      1 24   4 8.658914 38.60515 1.2682854   0.23
## 3   1      0 78   4 4.890728 28.22580 0.6015557   0.77
## 4   1      1 22   4 5.687345 42.62309 1.3919099   0.23
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 6.438753 0.5233818   4
## 2   1 100 5.073950 0.5630643   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##          method contrast  observed         se      lower       upper
## 1    asymptotic    A1-A0 -1.364803 0.76874573 -2.8715171  0.14191076
## 11    bootstrap    A1-A0 -1.364803 0.40418032 -1.5017040 -0.05010889
## 2    asymptotic    A1/A0  0.788033 0.10840007  0.6018039  1.03189094
## 2.5%  bootstrap    A1/A0  0.788033 0.09789247  0.6319720  0.98370057
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##          time        mcf     var_mcf     se_mcf arm
## 1 0.001067139 0.00000000 0.000000000 0.00000000   1
## 2 0.005045426 0.00000000 0.000000000 0.00000000   1
## 3 0.006134384 0.01045455 0.002295248 0.04790875   1
## 4 0.011786115 0.01045455 0.002295248 0.04790875   1
## 5 0.020344332 0.02090909 0.009180056 0.09581261   1
## 6 0.024756030 0.02090909 0.009180056 0.09581261   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##        method contrast  observed          p
## 1  asymptotic    A1-A0 -1.364803 0.07583787
## 3   bootstrap    A1-A0 -1.364803 0.03960396
## 5 permutation    A1-A0 -1.364803 0.07920792
## 2  asymptotic    A1/A0  0.788033 0.08331857
## 4   bootstrap    A1/A0  0.788033 0.03960396
## 6 permutation    A1/A0  0.788033 0.07920792
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate for continuous covariates. If covariates are provided, `CompareAUCs` uses an augmentation estimator to adjust for covariate differences.


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
  reps = 100,
  alpha = 0.05
)
show(aucs)
```

```
## Marginal Areas:
##   arm   n tau area    se
## 1   0 100   4 6.48 0.570
## 2   1 100   4 5.03 0.543
## 
## 
## CIs:
##       method contrast observed    se lower  upper
## 1 asymptotic    A1-A0    -1.49 0.775 -3.01 0.0274
## 2  bootstrap    A1-A0    -1.49 0.434 -1.63 0.0644
## 
## 
## P-values:
##        method contrast observed      p
## 1  asymptotic    A1-A0    -1.49 0.0543
## 2   bootstrap    A1-A0    -1.49 0.0990
## 3 permutation    A1-A0    -1.49 0.0990
```
