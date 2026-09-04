# Second-order variance implementation progress

## Current status

- Phase: 6. Documentation and final package check.
- Status: complete.
- Working tree was clean on branch `master` before implementation began.
- Generated Rcpp registration and roxygen documentation are synchronized with the completed implementation.

## Decisions and invariants

- `PseudoReg()` will be independent of the existing two-arm analysis functions.
- The public formula is one-sided and is evaluated on one row per subject.
- The regression is unweighted identity-link OLS fitted by `stats::lm()`.
- The corrected covariance uses denominator `n`, centered `h11 + h12` contributions, and normal Wald inference.
- MCF evaluation, risk sets, event ties, survival left limits, and jump weights must match `CalcMCF()` and `InfluenceFunction()`.
- The optimized projection will use design-weighted aggregate directions rather than allocating all subject pairs.
- Subjects in the derivative engine use the internal `FormatData()` order and are explicitly mapped back to the original identifiers before results are returned.
- The time grid is the ordered `CalcMCF()` grid through `tau`; risk indicators include a subject at its terminating time, survival is evaluated at the left limit, tied events are aggregated on that grid, and AUMCF increments have weight `(tau - u)`.

## Checkpoint record

### Checkpoint 1: shared pseudo-value preparation and subject alignment

- Completed by adding `.PreparePseudo()` and routing `GenPseudo()` through it without changing the public result.
- Original and internal identifiers are retained in `idx_map`; regression rows, pseudo-values, design rows, and derivative rows are explicitly matched by identifier.
- Existing pseudo-value tests passed unchanged.

### Checkpoint 2: reference R derivative implementation

- Completed in `tests/testthat/test-pseudo-regression.R` with a direct subject-pair implementation of the product-integral mixed second derivative.
- The reference covers MCF and AUMCF and serves as an independent oracle for the optimized implementation.

### Checkpoint 3: optimized Rcpp implementation

- Completed in `H12MCFCpp()` using design-weighted aggregate directions, with principal complexity `O(npm)` and no allocated `n` by `n` pair array.
- Rcpp registration is synchronized in `R/RcppExports.R` and `src/RcppExports.cpp`; the derivative engine remains internal.

### Checkpoint 4: public API and methods

- Completed with `PseudoReg()`, `print.PseudoReg()`, `summary.PseudoReg()`, `print.summary.PseudoReg()`, `coef.PseudoReg()`, `vcov.PseudoReg()`, and `confint.PseudoReg()`.
- Input validation covers the formula contract, confidence levels, time support, risk support, subject alignment, missing or changing covariates, jump weights, and full-rank model matrices.
- The returned `lm` fit retains conventional inference in `summary(object$fit)`, while all `PseudoReg` methods use corrected normal-Wald inference.

## Files changed

- `.Rbuildignore`: excludes the tracked `development` implementation journal from package builds.
- `development/second_order_variance_progress.md`: created this journal.
- `R/Pseudovalues.R`: added internal `.PreparePseudo()` and made `GenPseudo()` delegate to it while preserving its public return value.
- `src/PseudoRegression.cpp`: added the aggregate-direction empirical product-integral derivative engine.
- `R/RcppExports.R` and `src/RcppExports.cpp`: regenerated after adding `H12MCFCpp()`.
- `R/PseudoRegression.R`: added `PseudoReg()`, input/alignment validation, corrected covariance construction, and S3 methods.
- `tests/testthat/test-pseudo-regression.R`: added a direct pairwise R reference plus end-to-end, special-case, alignment, and failure-mode tests.
- `development/validate_second_order_variance.R`: added a fixed-seed validation driver for focused tests and repeated MCF/AUMCF centering checks.
- `README.Rmd`, `NEWS.md`, and `R/PseudoRegression.R`: added usage, release notes, assumptions, and an executable example.

## Checks run

- `git status --short`: clean before the two files above were added.
- Inspected `GenPseudo()`, `InfluenceFunction()`, `CalcMCF()`, `PsiMCF()`, `PsiAUC()`, package exports, and current pseudo-value tests.
- Existing pseudo-value tests: 19 passed, 0 failed; 39 pre-existing tidyselect warnings.
- `R CMD INSTALL --no-multiarch --with-keep.source . -l /tmp/mcc-lib`: derivative C++ compiled and package installed successfully.
- End-to-end smoke tests now produce finite corrected MCF and AUMCF covariance matrices; empirical `h12` column means were below `7e-18`.
- New test file first run: 14 assertions passed. The direct pairwise reference agreed for MCF and AUMCF. One formula-validation test exposed that `terms(~ .)` errors before the custom diagnostic; validation now detects `.` from the formula expression before calling `terms()`.
- New test file second run: 24 passed, 0 failed, 0 warnings.

## Known issues or open questions

- Initial nonfinite `h12` values were traced to and fixed by materializing an Armadillo expression before returning it from a local helper.
- The first finite-difference run showed that the existing `PsiMCF()` is the population-continuity IF evaluated with empirical quantities, not the exact jump-corrected derivative; the difference was finite-sample order and is permitted by the proof's score-negligibility condition. The derivative test now compares against the exact product-integral first derivative. The mixed-derivative identity-design extraction was corrected by its factor `n`. The exhausted-risk-set fixture now has two grid times so it reaches the intended derivative-domain check rather than an unrelated one-time-grid limitation in `CalcMCF()`.
- The distinction between the existing continuous-population empirical IF and the exact jump-corrected empirical derivative is documented here and covered by the finite-difference test; it does not leave an unresolved implementation issue.

## Exact next action

No implementation action remains. Future work is limited to the separately deferred simulation studies and possible post-fit adapters or model extensions.

## Resume entry: final validation

- Current phase and status: phase 5 is complete; phase 6 final package validation is in progress.
- Files and symbols changed in this step: `PseudoReg()` now enforces the lower as well as upper observed-time bound for `tau`; `confint.PseudoReg()` validates its requested confidence level; the failure-mode tests cover the lower bound.
- Commands and outcomes: the source-tarball `R CMD check --no-manual` completed with 308 passed assertions, no failures, and three warnings. One warning is an external R-header/compiler diagnostic and two concern the package's pre-existing vignette layout.
- Generated-file status: Rcpp registration and roxygen documentation are synchronized with the implemented public API; the final code-only validation above does not require regeneration.
- Known failures or unresolved questions: none in the estimator implementation. The final post-edit focused tests and source-package check remain to be run.
- Exact next action: reinstall the package, rerun the focused and fixed-seed validation tests, rebuild the source tarball, and repeat `R CMD check --no-manual`.

## Checkpoint 5: unit and numerical validation complete

- Current phase and status: phase 5 is complete; phase 6 documentation and package checking is the only remaining phase.
- Files and symbols changed: the special-case test now verifies the intercept-only and no-censoring conclusions independently; generated README output had two trailing spaces removed.
- Mathematical and API decisions confirmed: the optimized aggregate direction equals the direct pairwise average for MCF and AUMCF; empirical `h12` is centered; HC0 normalization is recovered when `h12` vanishes; the first and mixed second derivatives agree with finite differences.
- Commands and outcomes: a fresh package installation succeeded; `development/validate_second_order_variance.R` passed all 25 focused assertions and 50 fixed-seed regression fits. The maximum absolute `h12` column mean was `8.326673e-17`.
- Known failures or unresolved questions: none.
- Generated-file status: Rcpp registration and roxygen documentation are synchronized. README source and rendered Markdown are synchronized except for removal of output-only trailing whitespace.
- Exact next action: run `R CMD check --no-manual` on the freshly rebuilt `MCC_0.9.0.tar.gz`, inspect its final log, and record the final package status.

## Checkpoint 6: documentation and final package check complete

- Current phase and status: phase 6 and the implementation are complete.
- Files and symbols changed: `README.Rmd` and `README.md` contain a worked `PseudoReg()` example; `NEWS.md` records the feature; `man/PseudoReg.Rd` documents the API and assumptions; `.Rbuildignore` excludes `development` from package builds while retaining this journal in version control.
- Commands and outcomes: the final clean source tarball passed `R CMD check --no-manual` with 309 assertions, no failures, and three warnings. One warning is the host R-header/compiler mismatch for `-Wfixed-enum-extension`; two are the package's pre-existing absence of built `inst/doc` vignette copies. Examples, compiled code, R code, documentation, S3 consistency, tests, and vignette rebuilding all passed.
- Numerical tolerances and outcomes: optimized versus pairwise `h12` uses tolerance `2e-11`; first and mixed finite differences use `2e-7` and `2e-6`; zero and centering checks use `1e-12`. The repeated fixed-seed validation had maximum absolute `h12` column mean `8.326673e-17`.
- Generated-file status: `RcppExports`, `NAMESPACE`, and `man/PseudoReg.Rd` are synchronized with the final source. The source tarball includes the public R file, internal C++ engine, documentation, and tests, and excludes `development` as intended.
- Known failures or unresolved questions: none in this implementation. The three package-check warnings are unrelated to the new estimator.
- Exact next action: none; implementation is ready for review.

## Independent correctness audit: 2026-09-04

- Current phase and status: post-implementation mathematical and numerical audit complete; one row-alignment defect affecting nonconstant jump weights was found and corrected.
- Files and symbols changed: `R/Format.R` now attaches `jump_weights` before `ConvertIdxToInt()` can reorder rows; `tests/testthat/test-pseudo-regression.R` now verifies invariance to independent reordering of event and subject-level covariate rows.
- Mathematical checks: the C++ finite-product derivative was traced term by term to the proof; the aggregate direction was verified to equal the direct subject-pair average by bilinearity; the sign and factor on `h12` were independently checked against exact jackknife pseudo-values in a finite-population experiment.
- Independent finite-difference checks: across 10 generated datasets for each of MCF and AUMCF, the complete functional agreed with the package estimand to `8.9e-16`, the optimized `h12` agreed with Richardson-extrapolated mixed finite differences to `3.8e-8`, symmetry error was below `3.9e-9`, and column centering error was below `1.4e-16`. With tied times and nonunit jump weights, the corresponding maximum derivative discrepancy was `1.7e-8` and centering error was below `5.2e-17`.
- Projection check: for exact jackknife scores at sample sizes 40, 80, and 160, the RMSE against the projection using `h11 + h12` decreased from 0.268 to 0.127 for MCF and from 0.146 to 0.074 for AUMCF. Omitting, negating, or doubling `h12` left substantially larger, persistent score discrepancies, confirming the coefficient `+1` and excluding a missing factor of two.
- Empirical-IF convention: the existing package IF is the continuous-population formula evaluated empirically rather than the exact derivative of the finite empirical product integral. In an independent 20-replicate study, its `sqrt(n)`-scaled weighted-score discrepancy from the exact empirical derivative decreased approximately as `n^-1.10` for MCF and `n^-1.06` for AUMCF over `n = 50, 100, 200, 400`, supporting the score-negligibility condition used in the proof.
- HC3 simulation check: the v1 helper forms the outer product of `h11/(1-leverage) + h12`. This is a consistent HC3-style finite-sample variant, not a distinct asymptotic correction. An alternative additive convention, `HC3 + (corrected HC0 - HC0)`, changed cell-mean variances by at most 0.46% of HC3 across the 162 completed cells and did not alter the substantive finding that the second-order correction is small.
- Commands and outcomes: fresh `R CMD INSTALL --no-multiarch --with-keep.source` succeeded; the focused pseudo-regression test file passed all 28 assertions; the variance-calibration validation script passed; all remaining package tests passed when run against the installed package except five tests that call unexported helpers directly and therefore require `devtools::load_all()` rather than an attached installed namespace. The pseudo-regression test context itself had no failures.
- Defect impact: the jump-weight alignment defect did not affect the v0 or v1 simulations because they used unit jump weights. Before correction, a shuffled weighted dataset changed coefficients and covariance; after correction, the maximum discrepancies were `3.4e-16` for coefficients, `5.4e-16` for `h12`, and `1.4e-16` for covariance.
- Generated-file status: no Rcpp, namespace, roxygen, or TeX regeneration is required by the `Format.R` and test-only changes.
- Known failures or unresolved questions: no material error remains in the second-order estimator. The HC3-plus-second-order construction should be described explicitly when results are shared because HC3 is a finite-sample modification and several asymptotically equivalent extensions are possible.
- Final package check: a fresh source tarball passed `R CMD check --no-manual`, including all package tests and vignette rebuilding. Status was one warning caused by the host compiler not recognizing an R-header diagnostic flag; there were no estimator, test, code, or documentation failures.
- Exact next action: none; the estimator and the simulation conclusion are ready for substantive review.
