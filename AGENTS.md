# AGENTS.md

## What this package is

ertte is the time-to-event (TTE) member of the exposure-response (E-R)
package family, alongside [erglm](https://github.com/djnavarro/erglm)
(GLM E-R models) and [emaxnls](https://github.com/djnavarro/emaxnls)
(Emax / logistic-Emax models), all of which plug into
[erplots](https://github.com/djnavarro/erplots) for visualisation. See
[GitHub issue #1](https://github.com/djnavarro/ertte/issues/1) for the
full design/scoping discussion this package implements.

This first pass covers the issue's suggested phasing step 1 only:
**core modelling, methods, SCM, and simulation -- no plotting code**.
Two later phases from the issue are explicitly deferred (see "Planned
work" below):

- Phase 2: scalar E-R views of a TTE endpoint (landmark-binary / RMST)
  reusing erplots' existing `er_plot()`/`er_vpc()` grammars.
- Phase 3: a new `er_tte()` KM/survival-curve plotting grammar in the
  separate `erplots` repo.

`ertte_model()` wraps `survival::survreg()` -- parametric accelerated
failure time (AFT) model fitting (`ertte_model()`), survival-probability
prediction with confidence intervals (`ertte_predict()`), AIC-based
distribution selection (`ertte_select_distribution()`), stepwise
covariate modelling (`ertte_scm_forward()`/`ertte_scm_backward()`/
`ertte_scm_history()`, built on the single-term
`ertte_add_term()`/`ertte_remove_term()`), and simulation
(`ertte_fun()`, `simulate.ertte_model()`). Exponential, Weibull,
log-normal, and log-logistic distributions are tested and officially
supported.

The package's design is deliberately harmonised with the companion
`erglm` package -- closer to `erglm` than to `emaxnls`, since a fitted
`survreg` object (like a fitted `glm` object) already carries a rich set
of base S3 methods to inherit, unlike `emaxnls`'s bespoke NLS-based
class. `ertte_model()`, `ertte_predict()`, `ertte_fun()`,
`ertte_add_term()`/`ertte_remove_term()`, and the SCM
forward/backward/history functions all mirror their `erglm_*`
counterparts closely (see `R/erglm-core.R`/`R/erglm-scm.R` in the
companion repo). Genuine differences remain where the model classes
differ -- e.g. `survreg` significance testing is always a single
likelihood-ratio Chi-squared test (no family-dependent `test = "auto"`
argument, since `survreg`'s LRT doesn't vary by distribution the way
`glm()`'s does across families), and `survreg` doesn't retain its
fitting `data` on the returned object the way `glm()` does, so
`ertte_model()` stores it explicitly (`mod$data <- data`).

It deliberately contains **no plotting code**. For a model-agnostic
mini-language to visualise exposure-response models (including those
fitted here), see the companion package
[erplots](https://github.com/djnavarro/erplots). ertte interoperates
with erplots by implementing the `er_predict()`/`er_simulate()`/
`er_summary()` generics erplots defines, registered lazily at load time
(see `R/er-methods.R`) -- ertte has no hard dependency on erplots or on
plotting packages in package code.

## Planned work (deferred from the design issue, not yet done)

- **`coxph()` semi-parametric engine.** The design issue mentions this
  as optional ("Fit wrappers over base `survival` (`survreg` for
  parametric AFT; optionally `coxph` for a semi-parametric option)").
  Not implemented -- `ertte_model()` is `survreg`-only for now. A
  `coxph`-based model would need its own prediction/simulation story
  (no AFT location-scale structure to lean on) and probably a separate
  constructor rather than an extra `dist` value.
- **Power-function covariate parameterisation.** The design issue calls
  for "continuous covariates as power functions, categorical covariates
  as factors". `ertte_add_term()`/`ertte_remove_term()`/SCM currently
  only support plain additive terms on the AFT location scale (linear
  for continuous, factor levels for categorical) -- the erglm-style
  approach. The richer power-function parameterisation (closer to
  `emaxnls`'s per-structural-parameter covariate attachment, with
  delta-method or profile-likelihood CIs on the transformed parameter)
  is unimplemented.
- **`er_predict.ertte_model()`'s real contract.** Currently a minimal
  placeholder that forwards to `ertte_predict()` with a `time` argument
  threaded through `...` (not part of erplots' `er_predict(model,
  newdata, conf_level)` contract as currently defined, which assumes a
  response-vs-exposure view). Needs a proper design once phase 2
  (landmark-binary / RMST scalar E-R views, reusing `er_plot()`/
  `er_vpc()`) is scoped -- see the design issue's Workstream B.
- **Phase 3: `er_tte()` plotting grammar** -- lives in the separate
  `erplots` repo, co-designed with ertte per the issue. Not started.
- **Broader edge-case test coverage** from the issue's test-suite
  wishlist: all-censored data, single stratum, heavy ties, and the
  zero-exposure/placebo group specifically (though `ertte_data` does
  include a placebo group via `dose == 0`/`aucss == 0`).
- **Administrative-censoring simulation is a simplification.**
  `.ertte_simulate_draws()` caps every simulated event time at that
  row's *observed* exit time (`time`, whether that row was itself an
  event or a censoring), since the true administrative censoring time
  for subjects who had an event isn't otherwise available in a typical
  data set. A more accurate simulation would use a genuine per-subject
  administrative censoring/follow-up time where available, separately
  from the event indicator.

## Structure

- `R/ertte-core.R` -- `ertte_model()`, `ertte_predict()`, the
  `ertte_fun()` closure factory, and the shared
  `.ertte_simulate_draws()` helper (used directly by
  `er_simulate.ertte_model()` and by `simulate.ertte_model()` via
  `.ertte_resample()`). All four supported distributions are
  log-location-scale AFT models (`log(T) = mu + scale * W`); see
  `.ertte_dist_info()` in `R/utils-helpers.R` for the base
  distribution's CDF/quantile function this relies on.
- `R/ertte-family.R` -- `ertte_select_distribution()`: fits each
  candidate AFT distribution and returns the AIC-ranked comparison plus
  the best-fitting model.
- `R/ertte-scm.R` -- forward/backward stepwise covariate modelling
  (`ertte_scm_forward()`/`ertte_scm_backward()`/`ertte_scm_history()`),
  and the single-term `ertte_add_term()`/`ertte_remove_term()` helpers
  they're built on (also exported, matching erglm's
  `erglm_add_term()`/`erglm_remove_term()`).
- `R/ertte-simulate.R` -- `simulate.ertte_model()`, the `stats::simulate()`
  S3 method (and its `.ertte_resample()` helper), modelled on
  `simulate.erglm_model()`'s output shape: one row per observation per
  replicate, with `dat_id`/`sim_id`, sampled `coef_*` columns, and
  simulated `sim_time`/`sim_event`.
- `R/ertte-data.R` -- the synthetic `ertte_data` example dataset:
  structurally similar covariates to `erglm_data` (`sex`, `age`,
  `weight`, `dose`, `aucss`, `cmaxss`) for ecosystem consistency, plus
  simulated `time`/`event` from a Weibull AFT ground truth with an
  exposure and sex effect, administrative censoring at 180 days, and a
  little independent dropout censoring.
- `R/er-methods.R` -- erplots interoperability: S3 methods for
  `er_predict()`/`er_simulate()`/`er_summary()`, plus lazy registration
  via `.onLoad()` (vendored `.s3_register()`, copied from erglm's
  equivalent -- the standard pattern for optional cross-package S3
  methods).
- `R/utils-helpers.R`, `R/utils-global.R` -- small internal helpers,
  `globalVariables()` declarations for NSE, `.as_ertte()`, the
  `.ertte_check_*()` input validators, `.ertte_dist_info()` (the
  log-location-scale base-distribution table), and
  `.ertte_response_vars()` (extracts `time`/`event` variable names from
  a model's `Surv()` formula).

## Development workflow

- Document with roxygen2 (`devtools::document()`); Markdown roxygen is
  enabled (`Roxygen: list(markdown = TRUE)`).
- Run tests with `devtools::test()`; full checks with `devtools::check()`.
  The package currently checks cleanly (0 errors/warnings/notes).
- Tests live in `tests/testthat/`, one file per `R/` source file.
  `tests/testthat/test-er-methods.R` exercises interop with erplots and
  skips its erplots-registration check if erplots isn't installed.
- `ertte_data`'s generating code (`.make_ertte_data()`) is in
  `R/ertte-data.R`, commented out below the function definition (same
  pattern as `erglm_data`). Effect sizes/seed (`seed = 111L`,
  `beta_aucss = -0.0007`, `sex_effect = 0.4`) were tuned so that (a)
  `ertte_select_distribution()` recovers the true Weibull ground truth
  with a clear AIC margin, and (b) `ertte_scm_forward()` reliably
  detects the true `sex` effect at the default threshold -- both are
  exercised by regression tests, so if the dataset is ever regenerated
  with a different seed/effect size, check these tests still pass.

## Conventions

- Use the base R pipe (`|>`), not the magrittr pipe.
- Follow the existing tidyverse-style conventions (dplyr/tibble/rlang)
  already used throughout, matching erglm/emaxnls.
- Public functions are prefixed `ertte_`; internal helpers are prefixed
  with `.ertte_` (or, for a couple of package-wide utilities like
  `.pick_seed()`, no prefix at all).
- Model objects are plain `survreg` objects with an extra `ertte_model`
  class (same name as the constructor function `ertte_model()`,
  matching the base-R idiom of `lm()`/class `"lm"`) and an internal
  `$ertte` list for package-specific metadata (fitted `type`/`dist`, SCM
  history) -- see `.as_ertte()`.
- Don't add plotting code here -- that belongs in erplots.
