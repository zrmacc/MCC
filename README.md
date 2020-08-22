# Description

This package provides functions for inference on the difference and ratio in the areas under mean cumulative count curves comparing two treatment arms. 

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
* `time` is the event time. 
* `status` is coded 1 for the recurrent event, and 0 otherwise.
* `arm` is the treatment arm, coded as 1 for treatment, 0 for reference. 

For analyses of other data sets, arm and status should have the same coding. 

## Difference and Ratio of AUCs

To find a confidence interval and p-vaue for the difference and ratio in areas under the mean cumulative count curve at time $\tau = 60$: 

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
##   Time     Arm0     Arm1 Contrast    Estimate           L           U
## 1   60 115.3972 64.91554    A1-A0 -50.4816819 -78.8839927 -21.4368074
## 2   60 115.3972 64.91554    A1/A0   0.5625399   0.3863428   0.7844707
##            P
## 1 0.00149925
## 2 0.00199900
```
