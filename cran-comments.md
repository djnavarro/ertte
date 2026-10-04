## Resubmission (0.1.1)

This is a resubmission addressing two comments from the CRAN
reviewer on the 0.1 submission:

* "If there are references describing the methods in your package,
  please add these in the description field of your DESCRIPTION file
  ...". There are no standalone references describing the package's
  methods beyond the general time-to-event/survival-analysis
  literature already cited in the vignettes; `Surv()`/`survreg()`/
  `coxph()` themselves are documented (with their own references) in
  the `survival` package. No DESCRIPTION change was made for this
  comment.
* "Please add `\value` to .Rd files regarding exported methods ...
  -> Missing Rd-tags: reexports.Rd: `\value`". `R/reexports.R` now has
  a `@return` tag documenting `Surv()`'s return value (an object of
  class `"Surv"`), and the `@description` text was corrected to
  actually link to `survival::Surv()`'s documentation (it previously
  said "follow the links below" without any link being rendered).

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

* local Ubuntu 24.04, R 4.6.1, `devtools::check(cran = TRUE)` -- 0
  errors, 0 warnings, 0 notes (the `New submission` note from the 0.1
  cycle no longer appears locally, though CRAN's own incoming checks
  may still raise it since the package hasn't been accepted yet)
* R-hub v2 (`rhub::rhub_check()`), re-run for this 0.1.1 resubmission
  at commit `f4570cb`: linux, windows, macos-arm64, and nosuggests, all
  R-devel -- all four `Status: OK`
  (<https://github.com/djnavarro/ertte/actions/runs/37242824713>).
  macos-arm64 previously failed before `R CMD check` itself ran, due
  to rhub macOS runner/mirror infrastructure trouble (unreachable
  `bin/macos/arm64` binary repos forcing source compilation, which
  then failed loading the compiled `mvtnorm` dependency); that's
  resolved on this re-run.
* win-builder, re-run for this 0.1.1 resubmission via
  `devtools::check_win_release()`/`check_win_devel()`: R-release and
  R-devel both `Status: 1 NOTE` (the routine `New submission` note), 0
  errors, 0 warnings
  (<https://win-builder.r-project.org/69JuvIc5YjyZ/00check.log>,
  <https://win-builder.r-project.org/LtvbZNf1OXkv/00check.log>).
* ertte's regular CI also runs `R CMD check` via
  `r-lib/actions/check-r-package@v2` on `macos-latest`/`windows-latest`
  (R release) and three Linux R versions (devel/release/oldrel-1) on
  every push, as an additional cross-platform signal alongside rhub/
  win-builder.

## R CMD check results

0 errors | 0 warnings | 0 notes locally; CRAN's incoming checks may
still raise the routine `New submission` note, since the package
hasn't previously been accepted.

## Downstream dependencies

There are no downstream dependencies for this package (new submission).
