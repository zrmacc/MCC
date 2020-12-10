# Compare Mean Cumulative Count Curves

Zachary McCaw <br>
Updated: 2020-12-09



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
##   idx      time status arm      covar strat true_event_rate
## 1   1 0.6227499      2   1 -0.3260365     0       0.7438644
## 2   2 2.4232557      1   1  0.5524619     1       1.1311992
## 3   2 2.9472829      2   1  0.5524619     1       1.1311992
## 4   3 0.2746966      1   1 -0.6749438     0       0.6881469
## 5   3 0.9209057      0   1 -0.6749438     0       0.6881469
## 6   4 1.9251304      2   1  0.2143595     1       1.0489954
```

The essential data are:

* `idx`, the subject index. 
* `time`, the observation time. 
* `status`, coded 0 for censoring, 1 for an event, 2 for death.
* `arm`, coded as 1 for treatment, 0 for reference. 

For analyses of other data sets, arm and status should have the same coding. Each subject should experience at most one of censoring or death. Subjects who experience neither censoring nor death are assumed to remain at risk throughout follow-up. 

The example data also include:

* `covar`, a standard normal covariate that increases the recurrent event rate.
* `strat`, a two-level stratification factor that increases the event rate.
* `true_event_rate`, the patient-specific recurrent event rate.

### AUCs

To compare the areas under the mean cumulative count curves up to time $\tau = 4$: 

```r
aucs <- CompareStratAUCs(
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
##   Arm   N Tau Area    SE
## 1   0 100   4 6.23 0.686
## 2   1 100   4 5.09 0.551
## 
## 
## CIs:
##      Method            Type Contrast Observed     SE  Lower    Upper
## 1     Asymp      Equitailed    A1-A0   -1.140 0.8800 -2.870 0.582000
## 3 Bootstrap      Equitailed    A1-A0   -1.140 0.4830 -1.810 0.000763
## 4 Bootstrap Highest-density    A1-A0   -1.140 0.4830 -1.720 0.015000
## 2     Asymp      Equitailed    A1/A0    0.817 0.1260  0.603 1.110000
## 5 Bootstrap      Equitailed    A1/A0    0.817 0.0954  0.634 1.000000
## 6 Bootstrap Highest-density    A1/A0    0.817 0.0954  0.641 1.000000
##   Alpha_Lower Alpha_Upper
## 1      0.0250      0.0250
## 3      0.0250      0.0250
## 4      0.0303      0.0197
## 2      0.0250      0.0250
## 5      0.0250      0.0250
## 6      0.0298      0.0202
## 
## 
## P-values:
##      Method Contrast Observed      P
## 1     Asymp    A1-A0   -1.140 0.1940
## 3 Bootstrap    A1-A0   -1.140 0.0792
## 5      Perm    A1-A0   -1.140 0.1580
## 2     Asymp    A1/A0    0.817 0.1890
## 4 Bootstrap    A1/A0    0.817 0.0792
## 6      Perm    A1/A0    0.817 0.1580
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

The output of `CompareStratAUCs` is an object with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

```r
aucs@StratumAreas
```

```
##   arm strata   n tau     area var_area   se_area stratum_weight
## 1   0      1 100   4 6.231144 47.12374 0.6864673              1
## 2   1      1 100   4 5.087900 30.33698 0.5507901              1
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   Arm   N Tau     Area        SE
## 1   0 100   4 6.231144 0.6864673
## 2   1 100   4 5.087900 0.5507901
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##      Method            Type Contrast   Observed        SE      Lower
## 1     Asymp      Equitailed    A1-A0 -1.1432439 0.8801177 -2.8682429
## 3 Bootstrap      Equitailed    A1-A0 -1.1432439 0.4832584 -1.8075964
## 4 Bootstrap Highest-density    A1-A0 -1.1432439 0.4832584 -1.7201362
## 2     Asymp      Equitailed    A1/A0  0.8165274 0.1261156  0.6032532
## 5 Bootstrap      Equitailed    A1/A0  0.8165274 0.0953840  0.6335705
## 6 Bootstrap Highest-density    A1/A0  0.8165274 0.0953840  0.6413695
##          Upper Alpha_Lower Alpha_Upper
## 1 0.5817550597     0.02500     0.02500
## 3 0.0007625306     0.02500     0.02500
## 4 0.0149642365     0.03028     0.01972
## 2 1.1052027376     0.02500     0.02500
## 5 1.0001237979     0.02500     0.02500
## 6 1.0020909306     0.02980     0.02020
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##          time  mcf    var_mcf     se_mcf arm
## 1 0.001628056 0.01 0.00990000 0.09949874   1
## 2 0.011964590 0.02 0.01960000 0.14000000   1
## 3 0.019609591 0.03 0.02910000 0.17058722   1
## 4 0.047584250 0.03 0.02910000 0.17058722   1
## 5 0.060451115 0.04 0.03839992 0.19595897   1
## 6 0.075976616 0.05 0.06749980 0.25980723   1
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##      Method Contrast   Observed          P
## 1     Asymp    A1-A0 -1.1432439 0.19395522
## 3 Bootstrap    A1-A0 -1.1432439 0.07920792
## 5      Perm    A1-A0 -1.1432439 0.15841584
## 2     Asymp    A1/A0  0.8165274 0.18940745
## 4 Bootstrap    A1/A0  0.8165274 0.07920792
## 6      Perm    A1/A0  0.8165274 0.15841584
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

### Adjusted AUCs

The previous estimator allows for stratification, but a different approach is needed to accommodate for continuous covariates. The function `CompareAugAUCs` uses an augmentation estimator to adjusted for differences in 1 or more baseline covariates. 


```r
aucs <- CompareAugAUCs(
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
##   Arm   N Tau Area    SE
## 1   0 100   4 6.23 0.686
## 2   1 100   4 5.09 0.551
## 
## 
## CIs:
##      Method            Type Contrast Observed    SE Lower   Upper Alpha_Lower
## 1     Asymp      Equitailed    A1-A0    -1.15 0.877 -2.86  0.5720      0.0250
## 2 Bootstrap      Equitailed    A1-A0    -1.15 0.492 -1.84 -0.0882      0.0250
## 3 Bootstrap Highest-density    A1-A0    -1.15 0.492 -1.67  0.0195      0.0399
##   Alpha_Upper
## 1      0.0250
## 2      0.0250
## 3      0.0101
## 
## 
## P-values:
##      Method Contrast Observed      P
## 1     Asymp    A1-A0    -1.15 0.1910
## 2 Bootstrap    A1-A0    -1.15 0.0594
## 3      Perm    A1-A0    -1.15 0.1780
```
