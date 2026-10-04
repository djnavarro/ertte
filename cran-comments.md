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

* local Ubuntu 24.04, R 4.6.1, `devtools::check(cran = TRUE)`
* R-hub v2 (`rhub::rhub_check()`): linux, windows, and nosuggests, all
  R-devel -- all `Status: OK`. macos-arm64 (R-devel) failed before
  `R CMD check` itself ran: `setup-deps` reported the runner's
  macOS-arm64 binary CRAN/Bioconductor mirrors as unreachable
  (`pak::repo_status()` showed `ok = FALSE` for every `bin/macos/arm64`
  repo entry), which forced dependency compilation from source and
  then failed loading the compiled `mvtnorm` dependency (`dyn.load`).
  This is rhub macOS runner/mirror infrastructure trouble, not an
  issue with ertte -- corroborated by the two CRAN-infrastructure
  win-builder checks below both passing cleanly on the same code.
* win-builder, R-release and R-devel -- both `Status: 1 NOTE` (the
  routine `New submission` note), 0 errors, 0 warnings
  (<https://win-builder.r-project.org/9U18L0UUn2ND/00check.log>,
  <https://win-builder.r-project.org/57zjIk601rjZ/00check.log>)
* Since both rhub and CRAN's own macOS builder (`mac.r-project.org`,
  currently returning HTTP 502) are unavailable for this package right
  now, macOS coverage instead comes from ertte's regular CI, which
  runs `R CMD check` via `r-lib/actions/check-r-package@v2` on
  `macos-latest` (R release) on every push. The run at this exact
  commit (<https://github.com/djnavarro/ertte/actions/runs/35847857269>)
  passed cleanly on macOS, alongside Windows and three Linux R
  versions (devel/release/oldrel-1).

## R CMD check results

0 errors | 0 warnings | 1 note (`New submission`, expected for a
first submission)

## Downstream dependencies

There are no downstream dependencies for this package (new submission).
