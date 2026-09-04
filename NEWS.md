## Version 0.9.0

* Added `PseudoReg` for identity-link MCF and AUMCF pseudo-value regression with second-order-corrected covariance estimation.
* Added `CompareMvAUCs` for multivariate AUC contrasts with asymptotic and bootstrap inference.
* Added multivariate augmentation and process-specific eligibility weighting in influence functions.
* Extended `GenData` to simulate multiple recurrent event types with shared frailty.
* Added `PlotMvMCFs` and unified multivariate plotting via `PlotMCFs`.
* **Breaking change:** renamed univariate `weights` to `jump_weights` (argument and formatted data column).
* **Breaking change:** renamed multivariate contrast `weights` to `process_weights` in `CompareMvAUCs`.


## Version 0.8.4

* Reviewed code and added tests using Cursor.


## Version 0.8.3

* Updated augmentation weight.


## Version 0.8.2

* Added pseudo-value calcuation.


## Version 0.7.0

* Added the option to supply weights. 
* Added function for single-arm AUC calculation.


## Version 0.6.3

* Added ability to plot one sample mean cumulative function.
* Added plotting vignette.
