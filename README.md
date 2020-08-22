
# Description

This package provides functions for inference on the difference and ratio in the areas under mean cumulative count (MCC) curves, comparing two treatment arms. The MCC curves are estimated using the method of [Ghosh and Lin (2000)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.00554.x), which allows for the occurrence of terminal events such as death. 

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
##   idx time status arm
## 1   1    3      1   0
## 2   1   15      1   0
## 3   1   46      1   0
## 4   1   51      1   0
## 5   1   53      1   0
## 6   2   29      1   0
```

Here: 

* `idx` is the subject index. 
* `time` is the observation time. 
* `status` is coded 0 for censoring, 1 for an event, 2 for death.
* `arm` is the treatment arm, coded as 1 for treatment, 0 for reference. 

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
  reps = 2000,
  alpha = 0.05
)
show(aucs)
```

```
##   Time     Arm0     Arm1 Contrast    Estimate           L          U
## 1   60 103.5674 67.53429    A1-A0 -36.0330694 -62.0591162 -8.7798990
## 2   60 103.5674 67.53429    A1/A0   0.6520808   0.4632125  0.9019682
##             P
## 1 0.010494753
## 2 0.006996502
```
