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

`ertte_aft()` wraps `survival::survreg()` -- parametric accelerated
failure time (AFT) model fitting (`ertte_aft()`), survival-probability
prediction with confidence intervals (`ertte_predict()`), AIC-based
distribution selection (`ertte_select_distribution()`), stepwise
covariate modelling (`ertte_scm_forward()`/`ertte_scm_backward()`/
`ertte_scm_history()`, built on the single-term
`ertte_add_term()`/`ertte_remove_term()`), and simulation
(`ertte_fun()`, `simulate.ertte_model()`). Exponential, Weibull,
log-normal, and log-logistic distributions are tested and officially
supported. `ertte_coxph()` wraps `survival::coxph()` as a
semi-parametric sibling engine -- the constructor is scaffolded, but
its prediction/simulation story isn't implemented yet (see "Planned
work" and "API naming: AFT vs Cox PH" below).

The package's design is deliberately harmonised with the companion
`erglm` package -- closer to `erglm` than to `emaxnls`, since a fitted
`survreg` object (like a fitted `glm` object) already carries a rich set
of base S3 methods to inherit, unlike `emaxnls`'s bespoke NLS-based
class. `ertte_aft()`, `ertte_predict()`, `ertte_fun()`,
`ertte_add_term()`/`ertte_remove_term()`, and the SCM
forward/backward/history functions all mirror their `erglm_*`
counterparts closely (see `R/erglm-core.R`/`R/erglm-scm.R` in the
companion repo). Genuine differences remain where the model classes
differ -- e.g. `survreg` significance testing is always a single
likelihood-ratio Chi-squared test (no family-dependent `test = "auto"`
argument, since `survreg`'s LRT doesn't vary by distribution the way
`glm()`'s does across families), and `survreg`/`coxph` don't retain
their fitting `data` on the returned object the way `glm()` does, so
`ertte_aft()`/`ertte_coxph()` store it explicitly (`mod$data <- data`).

It deliberately contains **no plotting code**. For a model-agnostic
mini-language to visualise exposure-response models (including those
fitted here), see the companion package
[erplots](https://github.com/djnavarro/erplots). ertte interoperates
with erplots by implementing the `er_predict()`/`er_simulate()`/
`er_summary()` generics erplots defines, registered lazily at load time
(see `R/er-methods.R`) -- ertte has no hard dependency on erplots or on
plotting packages in package code.

## Planned work (deferred from the design issue, not yet done)

- **`coxph()` semi-parametric engine -- prediction/simulation.** The
  design issue mentions this as optional ("Fit wrappers over base
  `survival` (`survreg` for parametric AFT; optionally `coxph` for a
  semi-parametric option)"). The naming/API split is decided (see "API
  naming: AFT vs Cox PH" below), and `ertte_coxph()` (in
  `R/ertte-coxph.R`) is now scaffolded: it fits via `survival::coxph()`
  and returns an object with class `c("ertte_coxph", "ertte_model",
  "coxph")`, mirroring `ertte_aft()`. Still missing: `ertte_coxph`
  methods for `ertte_predict()`/`ertte_fun()`/`simulate()`. Unlike the
  closed-form `S(t)` available for AFT models, these need a baseline
  hazard estimate (e.g. via `survival::survfit()`), and simulation will
  need an analogous baseline-hazard-based counterpart to
  `.ertte_simulate_draws()`.
- ~~`ertte_add_term()`/`ertte_remove_term()` are AFT-hardcoded~~ --
  **fixed.** They used to refit by calling `ertte_aft(...)` directly
  (not `stats::update()`, which doesn't work here because the fitted
  object's `$call` refers to the constructor's own local argument
  bindings, not anything visible in the caller's frame). They now
  refit via an internal `.ertte_refit()` S3 generic (in
  `R/ertte-scm.R`) that dispatches on `mod`'s class
  (`.ertte_refit.ertte_aft()`/`.ertte_refit.ertte_coxph()`) and calls
  the matching constructor, so `ertte_scm_forward()`/
  `ertte_scm_backward()` now work for `ertte_coxph` models too. Fixing
  this also surfaced a genuine engine difference:
  `stats::drop.terms()`'s output (a `terms` object) can be passed
  straight back into `survreg()` as `formula`, but `coxph()` rejects it
  (errors in `terms.formula()`/`ExtractVars`) -- `ertte_remove_term()`
  now converts it via `stats::formula()` first, which works for both
  engines.
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

## API naming: AFT vs Cox PH

`ertte_aft()` (wraps `survreg()`) and `ertte_coxph()` (wraps `coxph()`)
are separate, engine-specific constructors -- not one constructor with
a `dist`/`engine` value -- since Cox PH is structurally different
enough (semi-parametric, no location-scale structure, no `dist`
argument) that folding it in would be misleading. There's no
back-compat alias for the old `ertte_model()` name; the package was
early enough in development that a clean rename was preferred.

Naming/dispatch scheme, now implemented for the constructors and
`ertte_predict()`/`ertte_fun()`, with `ertte_coxph`-specific
prediction/simulation still pending (see "Planned work" above):

- **Constructors are engine-specific by name**: `ertte_aft(formula,
  data, dist = "weibull", ...)` and `ertte_coxph(formula, data, ...)`.
  `dist` stays AFT-only; there's no Cox PH equivalent (nothing to
  select).
- **Shared superclass, engine-specific subclass**, for dispatch: AFT
  fits get class `c("ertte_aft", "ertte_model", "survreg")` (via
  `.as_ertte_aft()`), Cox PH fits get `c("ertte_coxph", "ertte_model",
  "coxph")` (via `.as_ertte_coxph()`), both in `R/utils-helpers.R`.
- **Downstream generics keep single shared names**, with S3 methods
  per subclass where behaviour genuinely differs:
  - `ertte_predict()` -- a generic (`UseMethod()`); only
    `ertte_predict.ertte_aft()` exists so far (closed-form `S(t)`).
    `ertte_predict.ertte_coxph()` (baseline-hazard-based) isn't
    implemented -- calling `ertte_predict()` on an `ertte_coxph` object
    currently errors with "no applicable method".
  - `ertte_fun()` -- same: a generic with only `ertte_fun.ertte_aft()`
    implemented so far.
  - `simulate()` -- only `simulate.ertte_model()` exists (used for AFT
    fits via the shared superclass); an `ertte_coxph`-specific method
    isn't implemented.
  - `ertte_scm_forward()`/`ertte_scm_backward()`/`ertte_scm_history()`/
    `ertte_add_term()`/`ertte_remove_term()` -- work across both
    engines. `ertte_add_term()`/`ertte_remove_term()` themselves stay
    single functions (not generics), but refit via an internal
    `.ertte_refit()` S3 generic that dispatches on `mod`'s class to
    call the matching constructor -- see the "Planned work" bullet
    above for the history of why a plain `stats::update()` doesn't work
    here.
  - `ertte_select_distribution()` -- stays AFT-only, no Cox PH
    equivalent.

## Structure

- `R/ertte-core.R` -- `ertte_aft()`, the `ertte_predict()`/`ertte_fun()`
  generics plus their `ertte_predict.ertte_aft()`/`ertte_fun.ertte_aft()`
  methods, and the shared `.ertte_simulate_draws()` helper (used
  directly by `er_simulate.ertte_model()` and by
  `simulate.ertte_model()` via `.ertte_resample()`; AFT-specific for
  now). All four supported distributions are log-location-scale AFT
  models (`log(T) = mu + scale * W`); see `.ertte_dist_info()` in
  `R/utils-helpers.R` for the base distribution's CDF/quantile function
  this relies on.
- `R/ertte-coxph.R` -- `ertte_coxph()`, the semi-parametric sibling
  constructor (wraps `survival::coxph()`). No `ertte_predict()`/
  `ertte_fun()`/`simulate()` methods yet -- see "Planned work" above.
- `R/ertte-family.R` -- `ertte_select_distribution()`: fits each
  candidate AFT distribution and returns the AIC-ranked comparison plus
  the best-fitting model.
- `R/ertte-scm.R` -- forward/backward stepwise covariate modelling
  (`ertte_scm_forward()`/`ertte_scm_backward()`/`ertte_scm_history()`),
  and the single-term `ertte_add_term()`/`ertte_remove_term()` helpers
  they're built on (also exported, matching erglm's
  `erglm_add_term()`/`erglm_remove_term()`). Both refit via the internal
  `.ertte_refit()` S3 generic (with `ertte_aft`/`ertte_coxph` methods),
  which is how they work across both engines despite not being
  generics themselves.
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
  `globalVariables()` declarations for NSE, the engine-specific class
  constructors `.as_ertte_aft()`/`.as_ertte_coxph()`, the
  `.ertte_check_*()` input validators, `.ertte_dist_info()` (the
  log-location-scale base-distribution table, AFT-only), and
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
- CI: `.github/workflows/R-CMD-check.yaml` (standard r-lib matrix:
  macOS/Windows/Ubuntu release, Ubuntu devel/oldrel-1) and
  `.github/workflows/pkgdown.yaml` (builds and deploys the pkgdown site
  to the `gh-pages` branch), both derived from the `r-lib/actions`
  examples and matching erglm's workflow files. Both install
  `github::djnavarro/erplots` as an extra dependency alongside CRAN
  packages, since `erplots` is a `Suggests`-only GitHub-hosted package
  exercised by `tests/testthat/test-er-methods.R`.
  `.github/workflows/test-coverage.yaml` runs the test suite under
  `covr` and uploads to Codecov (`CODECOV_TOKEN` secret configured on
  the repo), matching erglm's workflow. No `rhub.yaml` yet (erglm has
  one); add it later if/when CRAN submission becomes relevant.
- pkgdown renders every `*.md` file at the package root (and in
  `.github/`) into its own `docs/*.html` page -- hard-coded in
  `pkgdown:::package_mds()` and not configurable via `_pkgdown.yml`, so
  `.Rbuildignore`-ing `AGENTS.md` (needed to keep it out of the built
  *package*) has no effect on the *pkgdown site*: unhandled, it'd get
  published as `docs/AGENTS.html` and indexed in
  `docs/search.json`/`docs/sitemap.xml`. `tools/pkgdown-postbuild.R`
  strips this page (and its search/sitemap entries) back out;
  `.github/workflows/pkgdown.yaml` runs it right after
  `build_site_github_pages()`. Run it manually after any local
  `pkgdown::build_site()` too -- same pattern as erglm (which also
  strips a `PLAN.md`-derived page; ertte has no `PLAN.md`, so only
  `AGENTS.md` needs stripping here). The site uses the shared
  `djnavarro/waeponwifestre` template and a custom domain
  (`ertte.djnavarro.net`) matching erglm/emaxnls -- the DNS/GitHub Pages
  custom-domain configuration itself is done outside the repo.

## Conventions

- Use the base R pipe (`|>`), not the magrittr pipe.
- Follow the existing tidyverse-style conventions (dplyr/tibble/rlang)
  already used throughout, matching erglm/emaxnls.
- Public functions are prefixed `ertte_`; internal helpers are prefixed
  with `.ertte_` (or, for a couple of package-wide utilities like
  `.pick_seed()`, no prefix at all).
- Model objects are plain `survreg`/`coxph` objects with two extra
  classes prepended: an engine-specific subclass (`"ertte_aft"` or
  `"ertte_coxph"`, matching the constructor name, cf. the base-R idiom
  of `lm()`/class `"lm"`) ahead of the shared `"ertte_model"`
  superclass -- see `.as_ertte_aft()`/`.as_ertte_coxph()` in
  `R/utils-helpers.R`. Both carry an internal `$ertte` list for
  package-specific metadata (fitted `type`/`dist`, SCM history).
- Don't add plotting code here -- that belongs in erplots.
