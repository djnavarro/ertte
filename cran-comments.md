## Submission

This is a new release. ertte is the time-to-event member of the
exposure-response package family alongside
[erglm](https://github.com/djnavarro/erglm) (GLM-based models) and
[emaxnls](https://github.com/djnavarro/emaxnls) (Emax/logistic-Emax
models), both already on CRAN. It provides estimation, prediction, and
simulation tools for exposure-response models of time-to-event
endpoints built on `survival::survreg()`/`survival::coxph()`.

`erplots` (the model-agnostic visualisation package this family
interoperates with) is a `Suggests`-only dependency and is available
on CRAN (published 2026-09-09), so no `Additional_repositories` entry
is needed. Every use is conditional: its S3 methods are registered
lazily at load time only if `erplots` is present, and the one test
file exercising them is skipped via
`testthat::skip_if_not_installed("erplots")`.

## Test environments

* local Ubuntu 24.04, R 4.6.1, `devtools::check(cran = TRUE)`

R-hub and win-builder checks are planned as the next step in this
submission and will be added here before the package is uploaded.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies for this package (new submission).
