# Description

This package provides functions for inference on the difference and ratio in areas under mean cumulative count (MCC) curves, comparing two treatment arms. The MCC curves are estimated using the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which allows for the occurrence of terminal events such as death. 

# Installation


```r
devtools::install_github(repo = 'zrmacc/MCC')
```

# Examples

Synthetic example data in the format expected by this package may be loaded via:


```r
library(MCC)
data(mcc_data)
head(mcc_data)
```

```
##   idx time status arm strata
## 1   1    3      1   0      1
## 2   1   15      1   0      3
## 3   1   46      1   0      1
## 4   1   51      1   0      2
## 5   1   53      1   0      2
## 6   2   29      1   0      1
```

Here: 

* `idx` is the subject index. 
* `time` is the observation time. 
* `status` is coded 0 for censoring, 1 for an event, 2 for death.
* `arm` is the treatment arm, coded as 1 for treatment, 0 for reference. 
* `strata` is a 3-level stratification factor. 

For analyses of other data sets, arm and status should have the same coding. Each subject should experience at most one of censoring or death. Subjects who experience neither censoring nor death are assumed to remain at risk throughout follow-up. 

## Difference and Ratio of AUCs

To compare the areas under the mean cumulative count curves up to time $\tau = 60$: 

```r
set.seed(100)
aucs <- CompareAUCs(
  time = mcc_data$time,
  status = mcc_data$status,
  arm = mcc_data$arm,
  idx = mcc_data$idx,
  tau = 60,
  strata = mcc_data$strata,
  reps = 100,
  alpha = 0.05
)
show(aucs)
```

```
## Areas:
##    N0    Area0  N1    Area1
## 1 184 17.61869 143 14.79571
## 
## 
## CIs:
##            Method Contrast   Observed      Lower       Upper Lower_alpha
## 1       Equi-tail    A1-A0 -2.8229848 -5.8755215  0.28357670     0.02500
## 2 Highest-density    A1-A0 -2.8229848 -6.0697516 -0.05520067     0.01012
## 3       Equi-tail    A1/A0  0.8397733  0.6901222  1.01650562     0.02500
## 4 Highest-density    A1/A0  0.8397733  0.7072868  1.03006533     0.03988
##   Upper_alpha
## 1     0.02500
## 2     0.03988
## 3     0.02500
## 4     0.01012
## 
## 
## P-values:
##   Method Sides Contrast    P
## 1   Boot     1    A1-A0 0.04
## 2   Boot     1    A1/A0 0.04
## 3   Perm     1    A1-A0 0.07
## 4   Perm     2    A1-A0 0.11
## 5   Perm     1    A1/A0 0.06
## 6   Perm     2    A1/A0 0.09
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* The stratification factor `strata` is optional, and may be omitted. No strata should be empty in either arm.
* `reps` is the number of bootstrap replicates. The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable. 
* `alpha` is 1 minus the desired confidence interval (CI) coverage. 

### Outputs

The output of `CompareAUCs` is an object of class `compAUCs` with these slots.

* `@Areas` containing the overall sample size and AUCs for each arm:


```r
aucs@Areas
```

```
##    N0    Area0  N1    Area1
## 1 184 17.61869 143 14.79571
```

* `@CIs` containing the observed difference and ratio with confidence intervals:


```r
aucs@CIs
```

```
##            Method Contrast   Observed      Lower       Upper Lower_alpha
## 1       Equi-tail    A1-A0 -2.8229848 -5.8755215  0.28357670     0.02500
## 2 Highest-density    A1-A0 -2.8229848 -6.0697516 -0.05520067     0.01012
## 3       Equi-tail    A1/A0  0.8397733  0.6901222  1.01650562     0.02500
## 4 Highest-density    A1/A0  0.8397733  0.7072868  1.03006533     0.03988
##   Upper_alpha
## 1     0.02500
## 2     0.03988
## 3     0.02500
## 4     0.01012
```

* `@Curves` containing the per arm mean cumulative count curve, averaged across strata if applicable.


```r
head(aucs@Curves)
```

```
##   Time        MCF Arm
## 1    1 0.00000000   0
## 2    2 0.01760531   0
## 3    3 0.07808519   0
## 4    5 0.18608331   0
## 5    6 0.22229281   0
## 6    7 0.25892352   0
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##   Method Sides Contrast    P
## 1   Boot     1    A1-A0 0.04
## 2   Boot     1    A1/A0 0.04
## 3   Perm     1    A1-A0 0.07
## 4   Perm     2    A1-A0 0.11
## 5   Perm     1    A1/A0 0.06
## 6   Perm     2    A1/A0 0.09
```

* `@Reps` containing the bootstrap and permutation test statistics.


```r
head(aucs@Reps)
```

```
##       boot_diff boot_ratio  perm_diff perm_ratio
## [1,] -3.2083067  0.8295516  1.5749908  1.1012635
## [2,] -1.4661234  0.9237945 -0.2833097  0.9828161
## [3,] -2.4171223  0.8646791  0.6115070  1.0380623
## [4,] -1.8367062  0.8846476 -0.7210144  0.9567278
## [5,] -0.9420859  0.9445555  0.1193121  1.0073297
## [6,] -4.3356346  0.7749201  1.2256497  1.0775559
```

* `@Weights` containing the per-stratum weights and AUCs.


```r
aucs@Weights
```

```
##   Stratum Stratum_weight N0    Area0 N1    Area1
## 1       1      0.3211009 57 51.27138 48 34.58627
## 2       2      0.3302752 61 52.15020 47 54.58371
## 3       3      0.3486239 66 54.98438 48 43.75429
```
