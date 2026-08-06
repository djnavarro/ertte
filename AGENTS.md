# AGENTS.md

## What this package is

ertte is the time-to-event (TTE) member of the exposure-response (E-R)
package family, alongside [erglm](https://github.com/djnavarro/erglm)
(GLM E-R models) and [emaxnls](https://github.com/djnavarro/emaxnls)
(Emax / logistic-Emax models), all of which plug into
[erplots](https://github.com/djnavarro/erplots) for visualisation. See
[GitHub issue #1](https://github.com/djnavarro/ertte/issues/1) for the
full design/scoping discussion this package implements (the only issue
still open -- everything else scoped there has shipped; see
`.agents/PLAN.md` for what, if anything, remains).

ertte provides two model-fitting engines (`ertte_aft()`/`ertte_coxph()`),
prediction/simulation/scalar-reduction methods for both, stepwise
covariate modelling, and interoperability with erplots'
`er_predict()`/`er_simulate()`/`er_summary()`/`er_predict_survival()`
generics. It deliberately contains **no plotting code** -- that's
erplots' job.

## Architecture reference (current state)

This section documents how the package works *today*. For the design
rationale behind these choices, rejected alternatives, and a record of
how the API got here, see [.agents/HISTORY.md](.agents/HISTORY.md).

### Two engines, one shared superclass

`ertte_aft(formula, data, dist = "weibull", ...)` wraps
`survival::survreg()` -- a parametric, log-location-scale AFT model
(`log(T) = mu + scale * W`). Exponential, Weibull, log-normal, and
log-logistic distributions are supported;
`ertte_aft_select_distribution()` fits each candidate and returns the
AIC-ranked comparison plus the best model.

`ertte_coxph(formula, data, ...)` wraps `survival::coxph()` (fit with
`model = TRUE`, needed so `survfit()`/`basehaz()` don't try to
reconstruct the model frame from `object$call$data`, which fails --
see "Refitting without `stats::update()`" below) -- a semi-parametric
proportional-hazards model, `h(t|x) = h0(t) * exp(x'beta)`.

Both constructors are engine-specific by name (no `dist`/`engine`
switch on one constructor), since Cox PH has no location-scale
structure or `dist` argument to select over. Fitted objects carry two
extra classes prepended to the base engine class: `c("ertte_aft",
"ertte_model", "survreg")` or `c("ertte_coxph", "ertte_model",
"coxph")` (via `.as_ertte_aft()`/`.as_ertte_coxph()` in
`R/utils-helpers.R`), plus an internal `$ertte` list for
package-specific metadata (fitted `type`/`dist`, SCM history). Both
engines store their fitting `data` explicitly (`mod$data <- data`),
since neither `survreg`/`coxph` retains it the way `glm()` does.

Downstream generics keep single shared names, dispatching per engine
where behaviour genuinely differs:

- **`ertte_predict(object, newdata, time, conf_level, ...)`** -- point
  survival-probability predictions + CI. `ertte_predict.ertte_aft()`
  is closed-form (base-distribution CDF); `ertte_predict.ertte_coxph()`
  goes through `survival::survfit()` (baseline hazard + linear
  predictor, log-transform CI, `extend = TRUE` to hold survival
  constant beyond the last observed follow-up).
- **`ertte_fun(object, param = NULL)`** -- returns a closure,
  `function(data = NULL, time, param = NULL)`, evaluating a
  (possibly counterfactual, via `param`) survival curve without
  refitting. `ertte_fun.ertte_aft()` evaluates the closed-form `S(t)`
  directly; `ertte_fun.ertte_coxph()` evaluates
  `S0(t)^exp((x - xbar)'param)` via `survival::basehaz(object, centered
  = TRUE)` and `object$means`. A custom `param` only varies the linear
  predictor -- the baseline hazard/means are always taken from the
  fitted object, never recomputed (that would need refitting the
  partial likelihood's risk sets).
- **`simulate.ertte_model()`** -- the one `stats::simulate()` method
  for both engines; per-engine event-time simulation is dispatched one
  level down via the internal `.ertte_simulate_draws()` S3 generic
  (`.ertte_simulate_draws.ertte_aft()` in `R/ertte-aft.R`,
  `.ertte_simulate_draws.ertte_coxph()` in `R/ertte-coxph.R`). Takes an
  optional `censor_time` (a genuine per-row administrative follow-up
  time); absent that, censored rows are capped at their observed exit
  time and event rows are left uncensored (an approximation, but a
  less biased one than capping every row at its observed exit time,
  which would leak the observed outcome into the simulation).
- **`ertte_landmark(object, newdata, landmark_time, conf_level)`** --
  a single, non-generic function (not `ertte_predict()`/`ertte_fun()`
  style dispatch) that reduces a TTE endpoint to a binary landmark
  response, `P(event by t*) = 1 - S(t*)`, by calling `ertte_predict()`
  at `time = landmark_time` and swapping the CI bounds (a monotonic
  decreasing transform). Needs no engine-specific logic of its own.
- **`ertte_rmst(object, newdata, tau, conf_level)`** -- a genuine
  generic (unlike `ertte_landmark()`), since restricted mean survival
  time needs the whole curve integrated, not one `ertte_predict()`
  call. `ertte_rmst.ertte_aft()` integrates the closed-form `S(t|x)`
  via `stats::integrate()` (reparameterised onto the `log(t)` scale --
  see `.agents/HISTORY.md` for why) with an analytic delta-method SE
  (propagating only `Var(mu)`, matching `ertte_predict.ertte_aft()`'s
  own simplification). `ertte_rmst.ertte_coxph()` computes an *exact*
  step-function sum (not quadrature) for the point estimate, with an
  SE built on `survfit()`'s profile-specific `std.err` field (see
  `.agents/HISTORY.md` for why this isn't `survival:::survmean()`'s
  own SE, which understates uncertainty for leveraged covariate
  profiles). Both engines' `ci_lower`/`ci_upper` are symmetric Wald
  intervals on the RMST scale, not automatically bounded to `[0,
  tau]` -- see `.agents/PLAN.md`.
- **`ertte_scm_forward()`/`ertte_scm_backward()`/`ertte_scm_history()`**,
  built on the single-term `ertte_add_term()`/`ertte_remove_term()` --
  work across both engines via an internal `.ertte_refit()` S3 generic
  (dispatching on the model's class) rather than `ertte_add_term()`/
  `ertte_remove_term()` themselves being generics. A `criterion`
  argument (`"p-value"` (default, LRT-based, uses `threshold`), `"aic"`,
  or `"bic"`) selects the significance test; `ertte_scm_history()`'s
  history tibble records which criterion drove each step.
- **`ertte_power(x, ref = NULL)`** -- a power-function covariate
  transform, `log(x / ref)` (`ref` defaults to `median(x, na.rm =
  TRUE)`), usable as a model term (e.g. `ertte_power(age)`) in either
  engine's formula or an SCM candidate set. The fitted coefficient
  *is* the power exponent (`T = T_ref * (x/ref)^theta` for AFT,
  `h(t|x) = h0(t) * (x/ref)^theta` for Cox), and its ordinary Wald CI
  already is the CI on that exponent -- no delta method or profile
  likelihood needed. Every non-missing `x` must be strictly positive.
  A `makepredictcall.ertte_power()` method (the same mechanism
  `stats::poly()` uses) keeps predict-time evaluation consistent with
  the original fit's `ref`.

### Refitting without `stats::update()`

`survreg`/`coxph` fits carry a `$call` that refers to the constructor's
own local `formula`/`data` bindings, invisible in a caller's frame, so
`stats::update()` doesn't work on them. Every internal refit (SCM,
`ertte_add_term()`/`ertte_remove_term()`) goes through the
`.ertte_refit()` S3 generic (`R/ertte-scm.R`) instead, which calls the
matching constructor directly.

### erplots interoperability (`R/er-methods.R`)

Registered lazily at load time (vendored `.s3_register()`, matching
erglm's pattern) so ertte has no hard dependency on erplots.

- **`er_predict.ertte_model()`** forwards to `ertte_landmark()` (via a
  `landmark_time` argument through `...`) or `ertte_rmst()` (via
  `tau`) -- mutually exclusive, exactly one required, each renaming
  its fitted-value column to erplots' expected `fit_resp`.
- **`er_simulate.ertte_model()`** likewise takes `landmark_time`/`tau`
  through `...` and returns landmark- or RMST-transformed
  `fit_resp`/`sim_resp` draws (via the internal
  `.ertte_simulate_scalar_resp()`), for `er_vpc()` parity. A replicate
  whose simulated outcome is ambiguous relative to
  `landmark_time`/`tau` (censored before it) becomes `NA` -- the same
  complete-case convention a landmark/RMST analysis applies to
  genuinely censored observed data. See `.agents/PLAN.md` for the
  rejected IPCW/pseudo-value alternative.
- **`er_predict_survival.ertte_model()`** powers erplots'
  `er_tte_add_model()` (a parametric `S(t)` overlay on a Kaplan-Meier
  plot): a thin pass-through to `ertte_predict()` for every strictly
  positive `time_grid` entry, special-casing `time_grid == 0` (`S(0) =
  1` always, no model evaluation needed).
- `er_predict()`/`er_summary()`/`er_simulate()` all declare `...`,
  which erplots splices caller-supplied `predict_args`/`summary_args`/
  `simulate_args` into -- e.g. `er_plot_add_model(mod, predict_args =
  list(landmark_time = 90))`.

## Structure

- `R/ertte-aft.R` -- `ertte_aft()`; the `ertte_predict()`/`ertte_fun()`
  generics plus their `ertte_predict.ertte_aft()`/`ertte_fun.ertte_aft()`
  methods; the `.ertte_simulate_draws()` S3 generic (with its
  `ertte_aft` method here; the `ertte_coxph` method is in
  `R/ertte-coxph.R`). See `.ertte_dist_info()` in `R/utils-helpers.R`
  for the base distribution CDF/quantile/density table.
- `R/ertte-coxph.R` -- `ertte_coxph()`; `ertte_predict.ertte_coxph()`;
  `ertte_fun.ertte_coxph()` (using the internal
  `.ertte_coxph_basehaz_at()` step-function helper);
  `.ertte_simulate_draws.ertte_coxph()` (event-time simulation by
  inverting the fitted baseline cumulative hazard, via
  `.ertte_coxph_invert_basehaz()`).
- `R/ertte-rmst.R` -- `ertte_rmst()` generic, with
  `ertte_rmst.ertte_aft()` and `ertte_rmst.ertte_coxph()` (the latter
  via the internal `.ertte_rmst_pfun_delta()` helper). Has a
  website-only article, `vignettes/articles/rmst.Rmd`.
- `R/ertte-landmark.R` -- `ertte_landmark()`. Has its own website-only
  article, `vignettes/articles/landmark.Rmd`.
- `R/ertte-power.R` -- `ertte_power()`, the power-function covariate
  transform, plus its `makepredictcall()` method.
- `R/ertte-family.R` -- `ertte_aft_select_distribution()`.
- `R/ertte-scm.R` -- `ertte_scm_forward()`/`ertte_scm_backward()`/
  `ertte_scm_history()`, `ertte_add_term()`/`ertte_remove_term()`, and
  the internal `.ertte_refit()` S3 generic. Has a website-only article,
  `vignettes/articles/scm.Rmd`.
- `R/ertte-simulate.R` -- `simulate.ertte_model()` (and its
  `.ertte_resample()` helper), modelled on `simulate.erglm_model()`'s
  output shape.
- `R/ertte-data.R` -- the synthetic `ertte_data` example dataset (`sex`,
  `age`, `weight`, `dose`, `aucss`, `cmaxss`, `admin_censor`), simulated
  from a Weibull AFT ground truth with an exposure and sex effect,
  administrative censoring at 180 days, and some independent dropout
  censoring. Generating code (`.make_ertte_data()`) lives commented out
  below the function definition. `seed = 111L`, `beta_aucss = -0.0007`,
  `sex_effect = 0.4` were tuned so `ertte_aft_select_distribution()`
  recovers the true Weibull ground truth with a clear AIC margin and
  `ertte_scm_forward()` reliably detects the true `sex` effect at the
  default threshold -- both exercised by regression tests, so check
  them if the dataset is ever regenerated with a different seed/effect
  size.
- `R/er-methods.R` -- erplots interoperability (see above), plus lazy
  registration via `.onLoad()`.
- `R/utils-helpers.R`, `R/utils-global.R` -- small internal helpers,
  `globalVariables()` declarations for NSE, the engine-specific class
  constructors `.as_ertte_aft()`/`.as_ertte_coxph()`, the
  `.ertte_check_*()` input validators, `.ertte_dist_info()`, and
  `.ertte_response_vars()` (extracts `time`/`event` variable names from
  a model's `Surv()` formula).

## Vignette structure

`vignettes/articles/` (pkgdown-only, excluded from the built package
via `.Rbuildignore`) holds four articles, in reading order:
`overview.Rmd` (an `ertte_aft()`/`ertte_coxph()` primer, the intended
starting point for a new user), `landmark.Rmd`, `rmst.Rmd`, and
`scm.Rmd`. Render locally with `rmarkdown::render()` (against a freshly
`devtools::install()`-ed copy of the package, not a stale installed
one) to check an article knits before pushing -- pkgdown's build step
doesn't run in this repo's CI, only on deploy.

## Conventions

- Use the base R pipe (`|>`), not the magrittr pipe.
- Follow the existing tidyverse-style conventions (dplyr/tibble/rlang)
  already used throughout, matching erglm/emaxnls.
- Public functions are prefixed `ertte_`; internal helpers are prefixed
  with `.ertte_` (or, for a couple of package-wide utilities like
  `.pick_seed()`, no prefix at all).
- No back-compat aliases for renamed functions -- the package is
  GitHub-only/pre-CRAN, so renames are done as straight renames across
  `R/`, `tests/`, and vignettes.
- Don't add plotting code here -- that belongs in erplots.

## Development workflow

- Document with roxygen2 (`devtools::document()`); Markdown roxygen is
  enabled (`Roxygen: list(markdown = TRUE)`).
- Run tests with `devtools::test()`; full checks with `devtools::check()`.
  The package currently checks cleanly (0 errors/warnings/notes).
- Tests live in `tests/testthat/`, one file per `R/` source file.
  `tests/testthat/test-er-methods.R` exercises interop with erplots and
  skips gracefully if erplots isn't installed, or if it predates a
  given feature (e.g. `er_predict_survival()`/`predict_args`).
- CI: `.github/workflows/R-CMD-check.yaml` (standard r-lib matrix),
  `.github/workflows/pkgdown.yaml` (builds/deploys the pkgdown site),
  and `.github/workflows/test-coverage.yaml` (covr -> Codecov). All
  three install `erplots` as an extra dependency from GitHub, since
  it's a `Suggests`-only GitHub-hosted package exercised by
  `tests/testthat/test-er-methods.R`. `R-CMD-check.yaml`/
  `test-coverage.yaml` are currently pinned to erplots'
  `feat/er-tte-core-scaffolding` branch rather than its default
  branch (the only place `er_tte()`/`er_predict_survival()` exist so
  far); `pkgdown.yaml` installs the default branch. See
  `.agents/PLAN.md` for reverting this once that branch merges
  upstream.
- pkgdown renders every `*.md` file at the package root into its own
  `docs/*.html` page -- hard-coded in `pkgdown:::package_mds()`, not
  configurable via `_pkgdown.yml`. `tools/pkgdown-postbuild.R` strips
  the resulting stray `AGENTS.html` page (and its search/sitemap
  entries) back out; `.github/workflows/pkgdown.yaml` runs it right
  after `build_site_github_pages()`. Run it manually after any local
  `pkgdown::build_site()` too. The site uses the shared
  `djnavarro/waeponwifestre` template and a custom domain
  (`ertte.djnavarro.net`), matching erglm/emaxnls.

## Keeping this documentation current

This file (`AGENTS.md`) should stay a lean, current-state reference --
if a change makes something above inaccurate, update it in place
rather than appending a note about the change.

Two companion files in `.agents/` (also excluded from the built
package via `.Rbuildignore`) carry the parts that don't belong here:

- **[.agents/HISTORY.md](.agents/HISTORY.md)** -- a condensed record of
  completed design decisions and their rationale (what was tried,
  rejected, and why), for context in future sessions. When you finish a
  piece of nontrivial design work, add an entry here rather than
  growing this file with "used to be X, now Y" narrative.
- **[.agents/PLAN.md](.agents/PLAN.md)** -- scoped-out future work and
  deferred/open items. When you finish something listed there, move its
  write-up into `HISTORY.md` and remove it from `PLAN.md` rather than
  marking it "done" in place.
