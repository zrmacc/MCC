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

For analyses of other data sets, arm and status should have the same coding. Each subject should experience at most one of censoring or death. Subjects who experience neither censoring nor death are assumed to remain at risk throughout follow-up 

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
  alpha = 0.05,
  ci_type = 'ETI',
  pval_type = 'perm'
)
show(aucs)
```

```
##   Time     Arm0     Arm1 Contrast   Estimate CI_type           L         U
## 1   60 44.23623 52.89457    A1-A0 -8.6583336     ETI -17.7354610 0.5196249
## 2   60 44.23623 52.89457    A1/A0  0.8363096     ETI   0.6897325 1.0100526
##   P_type         P
## 1   perm 0.0990099
## 2   perm 0.0990099
```

Arguments include: 

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* The stratification factor `strata` is optional, and may be omitted. Strata empty in one arm or the other are allowed; these strata receive weight zero when combining areas across arms.
* `reps` is the number of bootstrap replicates. The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable. 
* `alpha` is 1 minus the desired confidence interval (CI) coverage. CIs are obtained via the percentile method, with probability $\alpha / 2$ allocated to each tail. 
* `ci_type` is either `'ETI'` for an equi-tailed interval, or `'HDI'` for the highest density interval.
* `pval_type` is either `'perm'` or `'boot'`. For `'perm'`, treatment assignments are permuted on each iteration, and the p-value is the proportion of the *null* statistics that are as or more extreme than the *observed* statistics. For `'boot'`, the p-value is twice the proportion of bootstrap replicates on which the sign of the difference is areas is reversed. 

