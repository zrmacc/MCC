# Compare Mean Cumulative Count Curves

Zachary McCaw <br>
Updated: 2020-12-12



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
##   idx      time status arm      covar strata true_event_rate
## 1   1 0.6227499      2   1 -0.3260365      0       0.7438644
## 2   2 2.4232557      1   1  0.5524619      1       1.1311992
## 3   2 2.9472829      2   1  0.5524619      1       1.1311992
## 4   3 0.2746966      1   1 -0.6749438      0       0.6881469
## 5   3 0.9209057      0   1 -0.6749438      0       0.6881469
## 6   4 1.9251304      2   1  0.2143595      1       1.0489954
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
aucs <- CompareAUCs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  idx = data$idx,
  strata = data$strata,
  tau = 4,
  boot = TRUE,
  perm = TRUE,
  reps = 100,
  alpha = 0.05,
  seed = 100
)
show(aucs)
```

```
## Marginal Areas:
##   arm   n area    se tau
## 1   0 100 3.81 0.402   4
## 2   1 100 3.41 0.337   4
## 
## 
## CIs:
##       method            type contrast observed    se  lower upper alpha_lower
## 1 asymptotic      equitailed    A1-A0   -0.393 0.525 -1.420 0.635      0.0250
## 3  bootstrap      equitailed    A1-A0   -0.393 0.312 -0.874 0.219      0.0250
## 4  bootstrap highest-density    A1-A0   -0.393 0.312 -0.869 0.221      0.0298
## 2 asymptotic      equitailed    A1/A0    0.897 0.130  0.675 1.190      0.0250
## 5  bootstrap      equitailed    A1/A0    0.897 0.112  0.701 1.100      0.0250
## 6  bootstrap highest-density    A1/A0    0.897 0.112  0.703 1.100      0.0298
##   alpha_upper
## 1      0.0250
## 3      0.0250
## 4      0.0202
## 2      0.0250
## 5      0.0250
## 6      0.0202
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.393 0.454
## 3   bootstrap    A1-A0   -0.393 0.297
## 5 permutation    A1-A0   -0.393 0.535
## 2  asymptotic    A1/A0    0.897 0.451
## 4   bootstrap    A1/A0    0.897 0.297
## 6 permutation    A1/A0    0.897 0.535
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
##   arm strata  n tau     area  var_area   se_area stratum_weight
## 1   0      0 85   4 4.023612 16.501195 0.4406038          0.825
## 2   0      1 15   4 2.775884 14.534342 0.9843557          0.175
## 3   1      0 80   4 3.485717 11.898402 0.3856553          0.825
## 4   1      1 20   4 3.064176  8.115189 0.6369925          0.175
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   arm   n     area        se tau
## 1   0 100 3.805259 0.4022501   4
## 2   1 100 3.411947 0.3371287   4
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##       method            type contrast   observed        se      lower     upper
## 1 asymptotic      equitailed    A1-A0 -0.3933124 0.5248437 -1.4219872 0.6353624
## 3  bootstrap      equitailed    A1-A0 -0.3933124 0.3116947 -0.8735538 0.2185036
## 4  bootstrap highest-density    A1-A0 -0.3933124 0.3116947 -0.8693152 0.2209075
## 2 asymptotic      equitailed    A1/A0  0.8966398 0.1297419  0.6752288 1.1906526
## 5  bootstrap      equitailed    A1/A0  0.8966398 0.1124192  0.7009945 1.0959922
## 6  bootstrap highest-density    A1/A0  0.8966398 0.1124192  0.7033427 1.0974199
##   alpha_lower alpha_upper
## 1      0.0250      0.0250
## 3      0.0250      0.0250
## 4      0.0298      0.0202
## 2      0.0250      0.0250
## 5      0.0250      0.0250
## 6      0.0298      0.0202
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##          time        mcf     var_mcf    se_mcf arm
## 1 0.001628056 0.01031250 0.008401465 0.0916595   1
## 2 0.011964590 0.02062500 0.016590234 0.1288031   1
## 3 0.019609591 0.03093750 0.024566309 0.1567364   1
## 4 0.047584250 0.03093750 0.024566309 0.1567364   1
## 5 0.060451115 0.04112109 0.032140689 0.1792782   1
## 6 0.075976616 0.05130469 0.056318333 0.2373148   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##        method contrast   observed         p
## 1  asymptotic    A1-A0 -0.3933124 0.4536224
## 3   bootstrap    A1-A0 -0.3933124 0.2970297
## 5 permutation    A1-A0 -0.3933124 0.5346535
## 2  asymptotic    A1/A0  0.8966398 0.4508539
## 4   bootstrap    A1/A0  0.8966398 0.2970297
## 6 permutation    A1/A0  0.8966398 0.5346535
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate for continuous covariates. If covariates are provided, `CompareAUCs` uses an augmentation estimator to adjust for covariate differences.


```r
aucs <- CompareAUCs(
  time = data$time,
  status = data$status,
  arm = data$arm,
  idx = data$idx,
  covar = data$covar,
  tau = 4,
  boot = TRUE,
  perm = TRUE,
  reps = 100,
  alpha = 0.05,
  seed = 100
)
show(aucs)
```

```
## Marginal Areas:
##   arm   n tau area    se
## 1   0 100   4 3.82 0.404
## 2   1 100   4 3.40 0.333
## 
## 
## CIs:
##       method            type contrast observed    se lower upper alpha_lower
## 1 asymptotic      equitailed    A1-A0   -0.419 0.520 -1.44 0.601       0.025
## 2  bootstrap      equitailed    A1-A0   -0.419 0.313 -1.08 0.249       0.025
## 3  bootstrap highest-density    A1-A0   -0.419 0.313 -0.83 0.302       0.040
##   alpha_upper
## 1       0.025
## 2       0.025
## 3       0.010
## 
## 
## P-values:
##        method contrast observed     p
## 1  asymptotic    A1-A0   -0.419 0.421
## 2   bootstrap    A1-A0   -0.419 0.376
## 3 permutation    A1-A0   -0.419 0.495
```
