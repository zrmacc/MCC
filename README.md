
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
  pval_type = '2-sided'
)
show(aucs)
```

```
##   Time     Arm0     Arm1 Contrast   Estimate           L        U          P
## 1   60 44.23623 52.89457    A1-A0 -8.6583336 -17.7383583 1.759781 0.05940594
## 2   60 44.23623 52.89457    A1/A0  0.8363096   0.6920833 1.036865 0.06930693
```

Arguments include: 

* `tau` is the truncation time, or the time up to which the AUC is calculated. 
* The stratification factor `strata` is optional, and may be omitted. Strata empty in one arm or the other are allowed; these strata receive weight zero when combining areas across arms.
* `reps` is the number of bootstrap replicates. The bootstrap is grouped by `idx`, and stratified by `strata`, if applicable. 
* `alpha` is 1 minus the desired confidence interval (CI) coverage. CIs are obtained via the percentile method, with probability $\alpha / 2$ allocated to each tail. 
* `pval_type` is either `'2-sided'` or `'1-sided'`. For `'2-sided'`, on each resample, the per-patient treatment assignments are randomized, and null values for the test statistics area calculated. The final p-value represents the proportion of resamples on which the null test statistics were as or more extreme than observed. For `pval_type = '1-sided'`, the treatment assignments are not randomized; instead, the p-value is calculated as the proportion of resamples on which the sign of the bootstrapped difference in areas disagrees with the observed sign. 

