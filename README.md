# Compare Mean Cumulative Count Curves

Zachary McCaw <br>
Updated: 20-11-18


### Description

This package provides functions for inference on the difference and ratio in AUCs comparing two mean cumulative count (MCC) curves. The MCC curves are estimated using the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which allows for the occurrence of terminal events such as death. Also see:
* [CICs](https://github.com/zrmacc/CICs) for comparing cumulative incidence curves. 

## Installation


```r
devtools::install_github(repo = 'zrmacc/MCC')
```

## Examples

### Data

Synthetic example data in the format expected by this package may be loaded via:


```r
library(MCC)
data(mcc_data)
head(mcc_data)
```

```
##   idx       time status arm
## 1   1 0.37908421      1   1
## 2   1 0.88876029      2   1
## 3   2 1.25659316      0   1
## 4   3 0.09017335      1   1
## 5   3 0.50027480      1   1
## 6   3 1.13766337      2   1
```

Here: 

* `idx` is the subject index. 
* `time` is the observation time. 
* `status` is coded 0 for censoring, 1 for an event, 2 for death.
* `arm` is the treatment arm, coded as 1 for treatment, 0 for reference. 

For analyses of other data sets, arm and status should have the same coding. Each subject should experience at most one of censoring or death. Subjects who experience neither censoring nor death are assumed to remain at risk throughout follow-up. See `? mcc_data` for details of the data generating process.

### AUCs

To compare the areas under the mean cumulative count curves up to time $\tau = 4$: 

```r
aucs <- CompareAUCs(
  time = mcc_data$time,
  status = mcc_data$status,
  arm = mcc_data$arm,
  idx = mcc_data$idx,
  tau = 4,
  reps = 100,
  boot = TRUE,
  perm = TRUE,
  alpha = 0.05,
  seed = 100
)
show(aucs)
```

```
## Marginal Areas:
##   arm   n tau area se_area
## 1   0 100   4 7.95   0.813
## 2   1 100   4 5.80   0.602
## 
## 
## CIs:
##      Method            Type Contrast Observed  Lower  Upper Lower_alpha
## 1     Asymp      Equitailed    A1-A0    -2.15 -4.130 -0.164      0.0250
## 2     Asymp      Equitailed    A1/A0     0.73  0.549  0.971      0.0250
## 3 Bootstrap      Equitailed    A1-A0    -2.15 -3.180 -0.438      0.0250
## 4 Bootstrap Highest-density    A1-A0    -2.15 -3.220 -0.730      0.0100
## 5 Bootstrap      Equitailed    A1/A0     0.73  0.534  0.911      0.0250
## 6 Bootstrap Highest-density    A1/A0     0.73  0.532  0.863      0.0101
##   Upper_alpha
## 1      0.0250
## 2      0.0250
## 3      0.0250
## 4      0.0400
## 5      0.0250
## 6      0.0399
## 
## 
## P-values:
##      Method Sides Contrast Observed      P
## 1     Asymp     2    A1-A0    -2.15 0.0338
## 2     Asymp     2    A1/A0     0.73 0.0307
## 3 Bootstrap     1    A1-A0    -2.15 0.0000
## 4 Bootstrap     1    A1/A0     0.73 0.0000
## 5      Perm     1    A1-A0    -2.15 0.0100
## 6      Perm     2    A1-A0    -2.15 0.0100
## 7      Perm     1    A1/A0     0.73 0.0100
## 8      Perm     2    A1/A0     0.73 0.0100
```

Here:

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* `boot` indicates to construct bootstrap confidence intervals. 
* `perm` indicates to perform permutation tests for the difference and ratio of AUCs.
* `reps` is the number of simulation replicates. 
  - The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable. 
* `alpha` is 1 minus the desired coverage for confidence intervals. 

#### Outputs

The output of `CompareAUCs` is an object of class `compAUCs` with these slots.

* `@StratumAreas` containing the stratum-specific AUCs for each arm.

```r
aucs@StratumAreas
```

```
##   arm stratum   n tau     area var_area   se_area
## 1   1       1 100   4 5.802334 36.18870 0.6015704
## 2   0       1 100   4 7.948094 66.05757 0.8127581
```

* `@MargAreas` containing the AUCs for each arm, marginalized over any strata. 


```r
aucs@MargAreas
```

```
##   arm   n tau     area   se_area
## 1   0 100   4 7.948094 0.8127581
## 2   1 100   4 5.802334 0.6015704
```

* `@CIs` containing confindence intervals for the difference and ratio of AUCs.


```r
aucs@CIs
```

```
##      Method            Type Contrast   Observed      Lower      Upper
## 1     Asymp      Equitailed    A1-A0 -2.1457602 -4.1276150 -0.1639054
## 2     Asymp      Equitailed    A1/A0  0.7300283  0.5487647  0.9711656
## 3 Bootstrap      Equitailed    A1-A0 -2.1457602 -3.1801687 -0.4380015
## 4 Bootstrap Highest-density    A1-A0 -2.1457602 -3.2199011 -0.7299726
## 5 Bootstrap      Equitailed    A1/A0  0.7300283  0.5344165  0.9113646
## 6 Bootstrap Highest-density    A1/A0  0.7300283  0.5316845  0.8626085
##   Lower_alpha Upper_alpha
## 1     0.02500     0.02500
## 2     0.02500     0.02500
## 3     0.02500     0.02500
## 4     0.01000     0.04000
## 5     0.02500     0.02500
## 6     0.01009     0.03991
```

* `@MCF` containing the per arm mean cumulative count curve, averaged across strata.


```r
head(aucs@MCF)
```

```
##          time        mcf    var_mcf    se_mcf Arm
## 1 0.001342556 0.00000000 0.00000000 0.0000000   0
## 2 0.014901890 0.01010101 0.01009998 0.1004987   0
## 3 0.016177535 0.02020202 0.01999384 0.1413996   0
## 4 0.018204774 0.03030303 0.02968157 0.1722834   0
## 5 0.019051486 0.04040404 0.05956927 0.2440682   0
## 6 0.019795095 0.05050505 0.06884476 0.2623828   0
```

* `@Pvals` containing the bootstrap and permutation p-values.


```r
aucs@Pvals
```

```
##      Method Sides Contrast   Observed          P
## 1     Asymp     2    A1-A0 -2.1457602 0.03383279
## 2     Asymp     2    A1/A0  0.7300283 0.03070434
## 3 Bootstrap     1    A1-A0 -2.1457602 0.00000000
## 4 Bootstrap     1    A1/A0  0.7300283 0.00000000
## 5      Perm     1    A1-A0 -2.1457602 0.01000000
## 6      Perm     2    A1-A0 -2.1457602 0.01000000
## 7      Perm     1    A1/A0  0.7300283 0.01000000
## 8      Perm     2    A1/A0  0.7300283 0.01000000
```

* `@Reps` is a list containing the bootstrap and permutation test statistics.

* `@Weights` containing the per-stratum weights and AUCs.


```r
aucs@Weights
```

```
##   stratum weight   n  n0  n1
## 1       1      1 200 100 100
```
