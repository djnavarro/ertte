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

- **`coxph()` semi-parametric engine -- now fully implemented.** The
  design issue mentions this as optional ("Fit wrappers over base
  `survival` (`survreg` for parametric AFT; optionally `coxph` for a
  semi-parametric option)"). The naming/API split is decided (see "API
  naming: AFT vs Cox PH" below), and `ertte_coxph()` (in
  `R/ertte-coxph.R`) is scaffolded: it fits via `survival::coxph()`
  (with `model = TRUE` -- see below) and returns an object with class
  `c("ertte_coxph", "ertte_model", "coxph")`, mirroring `ertte_aft()`.
  `ertte_predict.ertte_coxph()` is implemented (also in
  `R/ertte-coxph.R`), via `survival::survfit(object, newdata,
  conf.int = conf_level)`: this computes each row's survival curve from
  the fitted baseline hazard (Breslow/Efron, matching `object$method`)
  and linear predictor, evaluates it at `time` via `summary(...,
  extend = TRUE)` (letting `time` exceed the last observed follow-up,
  holding survival constant beyond it), and takes confidence intervals
  from `survfit()`'s own log-transform CI rather than a hand-rolled
  Wald interval -- genuinely different from `ertte_predict.ertte_aft()`
  methodologically, which is expected given the different model
  structure, not a bug. Getting `survfit()` to work at all required
  fitting `ertte_coxph()` with `model = TRUE`: `survfit.coxph()`
  otherwise tries to reconstruct the model frame by re-evaluating
  `object$call$data`, which fails for the same reason `stats::update()`
  fails on these fits (see the `.ertte_refit()` bullet below) -- the
  captured call refers to `ertte_coxph()`'s own local `formula`/`data`
  bindings, not anything visible in `survfit()`'s caller's frame.
  Storing the model frame directly sidesteps that.
  `ertte_fun.ertte_coxph()` is also implemented (same file): it
  returns a function evaluating `S(t | x) = S0(t)^exp((x - xbar)'param)`
  via `survival::basehaz(object, centered = TRUE)` (held constant
  beyond the last observed time, via the internal
  `.ertte_coxph_basehaz_at()` helper, mirroring
  `ertte_predict.ertte_coxph()`'s extrapolation) and `object$means`
  (the covariate means `coxph()` centers on, needed because
  `basehaz()`'s baseline is relative to that, not to `x = 0`). As with
  the AFT method, a custom `param` only varies the linear predictor --
  the baseline hazard/means are always taken from the fitted `object`,
  never recomputed for a hypothetical `param` (a `param`-dependent
  baseline would need refitting the partial likelihood's risk sets).
  One genuine engine difference surfaced here: Cox models have no
  intercept (it's absorbed into the baseline hazard), but
  `stats::model.matrix()` on the model's `terms()` adds one anyway
  since `terms()` doesn't record its absence -- `ertte_fun.ertte_coxph()`
  drops that column explicitly so `ncol(mm)` matches
  `length(coef(object))`. Finally, `simulate()` now works for
  `ertte_coxph` fits too, via a new `.ertte_simulate_draws.ertte_coxph()`
  method (`.ertte_simulate_draws()` was turned into an S3 generic,
  dispatching on the model's class -- `simulate.ertte_model()` and
  `er_simulate.ertte_model()` needed **no changes** to pick this up,
  since both already call `.ertte_simulate_draws()` directly). Event
  times are simulated by inverse-CDF sampling on the *cumulative
  hazard* scale rather than sampling from a parametric quantile
  function: `S(t | x) = u` rearranges to `H0(t) = -log(u) / exp(lp)`,
  so a new `.ertte_coxph_invert_basehaz()` helper inverts the fitted
  (step-function) baseline cumulative hazard at that target value --
  `Inf` if the target exceeds every observed hazard value, which then
  gets capped at the row's observed exit time downstream via the same
  `pmin(sim_time_raw, obs_time)` censoring convention the AFT method
  uses. Verified by prototyping: simulating many draws at the baseline
  covariate profile and comparing the empirical proportion surviving
  past a given time to the fitted `S0(t)` reproduces it closely.
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
- ~~Power-function covariate parameterisation~~ -- **implemented.** The
  design issue calls for "continuous covariates as power functions,
  categorical covariates as factors", flagging CIs on the transformed
  parameter (delta method vs profile likelihood) as an open question.
  Investigated `emaxnls` for a precedent first: it turns out `emaxnls`
  does *not* implement power functions either -- its covariate model is
  plain additive/linear per structural parameter, fit via `nls()` -- so
  this was new design, not a port. The key insight that shaped it: unlike
  `emaxnls`'s genuinely nonlinear structural parameters, both `ertte_aft()`
  (log-location-scale AFT: `log(T) = mu + scale * W`) and `ertte_coxph()`
  (Cox PH: `h(t|x) = h0(t) * exp(x'beta)`) are already linear in their
  covariates on the model's natural (log-time / log-hazard-ratio) scale.
  A power-function effect -- `T = T_ref * (x/ref)^theta` (AFT) or `h(t|x)
  = h0(t) * (x/ref)^theta` (Cox) -- is, after taking logs, exactly a
  linear term in `log(x/ref)`. So `ertte_power(x, ref = NULL)` (in
  `R/ertte-power.R`, exported) just returns `log(x/ref)`: the fitted
  coefficient on `ertte_power(age)` *is* the power exponent directly, and
  its ordinary Wald CI from `confint()`/`summary()` already is the CI on
  that exponent -- no delta method or profile likelihood needed, which
  resolves the design issue's open question for these two model families
  specifically. `ref` defaults to `median(x, na.rm = TRUE)` in the fitting
  data (pop-PK/NONMEM convention); every non-missing `x` must be strictly
  positive (`log(x/ref)` is undefined otherwise), which rules out using it
  on covariates with a placebo/zero group (`dose`/`aucss`/`cmaxss`) -- by
  design, since the design issue's power-function language targets the
  *covariate model* (age, weight, ...), not the primary exposure metric.
  Predict-time consistency (reusing the original fit's `ref`, not
  recomputing a new median from `newdata`) is handled by a
  `makepredictcall.ertte_power()` method -- the same mechanism
  `stats::poly()`/`splines::ns()` use -- which required no changes to any
  of the model-matrix-building code elsewhere in the package
  (`ertte_predict()`/`ertte_fun()`/`.ertte_simulate_draws()`, either
  engine): they all already build design matrices via `model.matrix()`/
  `predict()`/`survfit()` against `object$terms` (or a `delete.response()`
  of it) on a plain `data.frame`, and `model.matrix.default()` always
  routes that through `model.frame()`, which honours a `terms` object's
  `predvars` attribute (populated at original-fit time from each term's
  `makepredictcall()`) -- confirmed empirically for both engines,
  including `survfit.coxph()`'s newdata path. Likewise,
  `ertte_add_term()`/`ertte_remove_term()`/SCM needed **no changes**:
  their term handling (`.ertte_check_term()`, `.ertte_check_candidates()`,
  `all.vars()`, `stats::drop.terms()`) already operates generically on
  formula term-labels rather than a restricted "bare variable name"
  grammar, so `~ ertte_power(age)` / `"ertte_power(age)"` were already
  legitimate terms/candidates once `ertte_power()` existed. Nothing
  prevents combining a plain linear term (`age`) and a power term
  (`ertte_power(age)`) for the same variable in the same model or
  candidate set -- deliberately left to the user's judgement, consistent
  with how term handling elsewhere in ertte doesn't reason about variable
  semantics.
- **`er_predict.ertte_model()`'s real contract -- landmark-binary now
  implemented; RMST and VPC parity still deferred.** Phase 2's
  Workstream B1 (scalar E-R views of a TTE endpoint) is partially done:
  a new exported `ertte_landmark(object, newdata, landmark_time,
  conf_level)` (in `R/ertte-landmark.R`) reduces a TTE endpoint to a
  binary landmark response, `P(event by t*) = 1 - S(t*)`, by calling
  `ertte_predict()` at `time = landmark_time` and transforming its
  survival-probability output -- since that's a decreasing monotonic
  transform, the confidence interval bounds simply swap (no
  recomputation, and no new edge-case handling needed: whatever
  validity `ertte_predict()`'s interval has for a given engine carries
  through unchanged). `ertte_landmark()` is deliberately **not** a
  generic, unlike `ertte_predict()`/`ertte_fun()` -- it needs no
  engine-specific logic of its own, since it delegates entirely to
  `ertte_predict()`, which already dispatches on the
  `ertte_aft`/`ertte_coxph` subclass. `er_predict.ertte_model()` (in
  `R/er-methods.R`) now forwards to it, threading a required
  `landmark_time` argument through `...` (erplots' `er_predict(model,
  newdata, conf_level)` contract has a fixed signature with no
  TTE-specific argument, so this -- like the old placeholder's `time`
  argument -- has to travel through `...`; omitting it errors
  informatively rather than silently doing nothing).

  Two pieces from the design issue's Workstream B1 are still
  deliberately deferred, confirmed with the maintainer when scoping
  this:
  - **RMST** (restricted mean survival time, the other scalar
    reduction the issue mentions) -- harder than landmark-binary, since
    it's an integral of `S(t)` over `[0, tau]`: the AFT engine has a
    closed-form `S(t)`, so an analytic RMST/CI is plausible, but the
    Cox engine's baseline hazard is an empirical step function, so its
    RMST would need numerical integration plus a delta-method or
    bootstrap CI -- no small addition.
  - **`er_simulate.ertte_model()` landmark-VPC parity.** The issue notes
    a landmark-binary VPC "likewise reuses `er_vpc()` unchanged", but
    that needs `er_simulate.ertte_model()` to also return
    landmark-transformed draws (`fit_resp`/`sim_resp` as event-by-t*
    indicators/probabilities), not its current per-row
    `sim_time`/`sim_event` shape (unchanged by this pass).

  **Caveat discovered end-to-end testing this against erplots'
  `er_plot()` (not an ertte bug -- tracked upstream as
  [erplots#10](https://github.com/djnavarro/erplots/issues/10)):**
  `er_plot_add_model(mod, landmark_time = 90)` currently errors, even
  though `er_predict.ertte_model()`'s contract (above) is implemented
  correctly and works when called directly. The reason lives entirely
  in `erplots`: `er_plot_add_model()`'s `...` is captured only for its
  *style builder* (`config$dots`, per `?er_style`), never forwarded to
  `er_predict()` itself -- `.get_model_predictions()` always calls
  `er_predict(model, newdata, conf_level)` with no extra arguments, so
  a required `landmark_time` can never reach it this way. The identical
  gap affects `er_plot_add_summary()`/`er_summary()` (which currently
  forwards *no* arguments at all, not even `conf_level`) and
  `er_vpc_add_simulated()`/`er_simulate()`. erplots#10 tracks the fix,
  with a design sketch comparing two options: forwarding the same
  `...` to both the style builder and the generic (risks silent
  argument collisions between the two), vs. a new, separate
  `predict_args`/`summary_args`/`simulate_args` argument per
  `er_plot_add_*()`/`er_vpc_add_*()` function (unambiguous, but more
  API surface). Until resolved upstream, plotting an `ertte_landmark()`
  curve via `er_plot_add_model()` needs a small wrapper object that
  bakes `landmark_time` into the model rather than passing it at
  `er_plot_add_model()` call time:

  ```r
  new_landmark_model <- function(object, landmark_time) {
    structure(list(object = object, landmark_time = landmark_time), class = "ertte_landmark_model")
  }
  er_predict.ertte_landmark_model <- function(model, newdata, conf_level = 0.95, ...) {
    ertte_landmark(model$object, newdata = newdata, landmark_time = model$landmark_time, conf_level = conf_level)
  }
  ```

  This workaround isn't currently shipped in the package (it's just a
  documented pattern here) -- revisit once erplots#10 lands, and
  reconsider whether ertte should export a constructor like this
  itself in the meantime.
- **Phase 3: `er_tte()` plotting grammar** -- lives in the separate
  `erplots` repo, co-designed with ertte per the issue. Not started.
- ~~Broader edge-case test coverage~~ -- **addressed** for the three
  named cases from the issue's test-suite wishlist (all-censored data,
  single stratum, heavy ties -- the zero-exposure/placebo group was
  already covered via `ertte_data`'s `dose == 0`/`aucss == 0` rows, used
  throughout the existing test suite). New tests live in
  `tests/testthat/test-ertte-edge-cases.R`. Findings, in order of how
  much they mattered:
  - **All-censored data (0 events) genuinely broke `ertte_coxph()`'s
    downstream methods.** `ertte_aft()` degrades gracefully (a
    `survreg()` "did not converge" warning, but `ertte_predict()` still
    returns sensible values). `ertte_coxph()` fit fine (`NA`
    coefficients, `n=0` events), but `ertte_predict()`/`ertte_fun()`/
    `simulate()` on it used to fail with a cryptic
    `model.frame.default()` error ("'data' must be a data.frame,
    environment, or list") -- traced to a `survival::coxph()` quirk: with
    zero events, `coxph(model = TRUE)` still leaves `object$model`
    unset, defeating the `model = TRUE` workaround `ertte_coxph()`
    otherwise relies on for `survfit()`/`basehaz()` to avoid
    reconstructing the model frame from `object$call$data` (which fails
    for the reason described below under "`ertte_add_term()`/
    `ertte_remove_term()`"). Fixed with a new
    `.ertte_check_coxph_nevent()` guard (in `R/ertte-coxph.R`), called at
    the top of `ertte_predict.ertte_coxph()`/`ertte_fun.ertte_coxph()`/
    `.ertte_simulate_draws.ertte_coxph()`, that aborts with an
    informative message instead. The constructor itself is left alone --
    fitting on all-censored data is still allowed (a legitimate, if
    degenerate, thing to inspect via `coef()`/`summary()`); only the
    baseline-hazard-based methods are actually unsupported.
  - **Single-level ("single stratum") categorical covariates** (e.g.
    `sex` with only one level present in a subset) fit fine on both
    engines, reporting the aliased level's coefficient as `NA`
    ("singularities" per `survreg()`'s own message). Downstream behaviour
    genuinely differs by *method*, not just by engine, and this is now
    pinned down by tests rather than left as a surprise: `predict.survreg()`
    (used by `ertte_predict.ertte_aft()`) propagates the `NA` straight
    through to every survival probability/CI; `survival::survfit()`
    (used by `ertte_predict.ertte_coxph()`) silently drops the aliased
    column and returns finite predictions anyway; `ertte_fun.ertte_coxph()`'s
    manual `mm %*% param` propagates the `NA` like the AFT method does.
    None of this crashes, and it isn't treated as a bug to fix -- SCM
    already handles the case where a candidate's addition produces a
    literal `NA` `anova()` p-value (skipped with a warning, pre-existing
    behaviour); a single-level factor candidate typically doesn't even
    trigger that path, since `anova()` reports a trivial `p = 1` (zero
    deviance) rather than `NA` when there's no variation left to explain.
  - **Heavy ties in event times** (many events at the same observed
    time) caused no issues on either engine -- `coxph()`'s default Efron
    tie-handling and `.ertte_simulate_draws()` both work unchanged. Kept
    as regression tests rather than a design concern.
- ~~Administrative-censoring simulation is a simplification~~ --
  **refined.** `.ertte_simulate_draws()` used to cap every simulated
  event time at that row's *observed* exit time (`time`), whether that
  row was itself an event or a censoring. Splitting by the observed
  `event` indicator sharpens the diagnosis: for censored rows, the
  observed exit time genuinely *is* the row's true censoring time (an
  exact match, not an approximation); for event rows, it isn't -- it's
  when the event happened, not the (necessarily later, unobserved)
  administrative censoring horizon -- so capping there leaked the
  observed outcome into the simulation and biased simulated-vs-observed
  comparisons (e.g. a VPC) toward looking more similar than the fitted
  model actually implies.
  - `simulate.ertte_model()`/`.ertte_simulate_draws()` (both engines) now
    take an optional `censor_time` argument (a single number, recycled,
    or a numeric vector of length `nrow(newdata)`) giving a genuine
    per-row administrative follow-up time, validated by a new
    `.ertte_check_censor_time()`. When supplied, it caps *every* row
    uniformly, regardless of observed event status -- the accurate case.
  - Absent `censor_time` (the default), the refined fallback -- via a
    new shared `.ertte_apply_admin_censoring()` helper -- caps censored
    rows at their observed exit time (unchanged) but leaves event rows
    **uncensored**, removing the specific bias above (still an
    approximation for event rows, just a less wrong one). This is a
    genuine default-behaviour change, confirmed with the maintainer
    before implementing.
  - `ertte_data` gained a demonstration column, `admin_censor` (fixed at
    180 for every row, reflecting the fixed study-wide cutoff
    `.make_ertte_data()` already used internally but didn't retain).
    Adding it as a plain constant (not an RNG draw) didn't perturb any
    other column's simulated values under the same `seed = 111L` --
    confirmed by regenerating and diffing before saving.
  - One real bug surfaced while implementing the coxph engine's version:
    `.ertte_coxph_invert_basehaz()` returns `Inf` whenever a simulated
    draw would need to survive past the fitted baseline hazard's
    support (the last observed follow-up time across the whole cohort)
    to "fail" -- previously always silently absorbed by the old blanket
    `pmin(sim_time_raw, obs_time)` cap, but left as a bare `Inf`
    (`sim_event = 1`, nonsensically "an event at time infinity") once
    event rows stopped being capped by default. Fixed by substituting
    `max(bh$time)` for `Inf` before applying `censor_time`/`obs_time`
    (so a smaller genuine cap still takes precedence), then forcing
    `sim_event = 0` for any row that hit this extrapolation boundary --
    consistent with the flat baseline-hazard extrapolation
    `ertte_predict.ertte_coxph()` already uses beyond the observed
    range. The AFT engine's equivalent (`qbase()` at `u` near 1) has the
    same theoretical failure mode but was left unguarded -- the risk is
    negligible there (a continuous quantile function, vs. the Cox
    engine's empirical/step-function baseline hazard, which realistically
    exhausts its support whenever `nsim` is reasonably large).

## API naming: AFT vs Cox PH

`ertte_aft()` (wraps `survreg()`) and `ertte_coxph()` (wraps `coxph()`)
are separate, engine-specific constructors -- not one constructor with
a `dist`/`engine` value -- since Cox PH is structurally different
enough (semi-parametric, no location-scale structure, no `dist`
argument) that folding it in would be misleading. There's no
back-compat alias for the old `ertte_model()` name; the package was
early enough in development that a clean rename was preferred.

Naming/dispatch scheme, now fully implemented for both engines
(constructors, `ertte_predict()`, `ertte_fun()`, and simulation -- see
"Planned work" above for how each landed):

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
  - `ertte_predict()` -- a generic (`UseMethod()`), now with methods
    for both engines: `ertte_predict.ertte_aft()` (closed-form `S(t)`)
    and `ertte_predict.ertte_coxph()` (baseline-hazard-based, via
    `survival::survfit()`).
  - `ertte_fun()` -- same: a generic, now with methods for both
    engines: `ertte_fun.ertte_aft()` and `ertte_fun.ertte_coxph()`
    (baseline-hazard-based, via `survival::basehaz()`).
  - `simulate()` -- only `simulate.ertte_model()` exists (used for both
    engines via the shared superclass) -- no separate
    `simulate.ertte_aft()`/`simulate.ertte_coxph()` methods, since the
    per-engine mechanics are dispatched one level down, by the internal
    `.ertte_simulate_draws()` S3 generic (see `R/ertte-core.R`).
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
  methods, and the `.ertte_simulate_draws()` S3 generic (with its
  `ertte_aft` method here; the `ertte_coxph` method is in
  `R/ertte-coxph.R`) -- used directly by `er_simulate.ertte_model()` and
  by `simulate.ertte_model()` via `.ertte_resample()`, both of which
  work for either engine automatically via this dispatch. All four
  supported distributions are log-location-scale AFT models (`log(T) =
  mu + scale * W`); see `.ertte_dist_info()` in `R/utils-helpers.R` for
  the base distribution's CDF/quantile function this relies on.
- `R/ertte-coxph.R` -- `ertte_coxph()`, the semi-parametric sibling
  constructor (wraps `survival::coxph()`, with `model = TRUE`),
  `ertte_predict.ertte_coxph()` (baseline-hazard-based survival
  predictions via `survival::survfit()`), `ertte_fun.ertte_coxph()`
  (counterfactual survival-curve evaluation via `survival::basehaz()`,
  using the internal `.ertte_coxph_basehaz_at()` step-function helper),
  and `.ertte_simulate_draws.ertte_coxph()` (event-time simulation by
  inverting the fitted baseline cumulative hazard, via
  `.ertte_coxph_invert_basehaz()`).
- `R/ertte-landmark.R` -- `ertte_landmark()`: the landmark-binary
  scalar E-R reduction (`P(event by t*) = 1 - S(t*)`) that
  `er_predict.ertte_model()` (`R/er-methods.R`) forwards to. A single
  function, not a generic -- it delegates entirely to
  `ertte_predict()`, which already dispatches per engine. See "Planned
  work" above for what's deferred (RMST, `er_simulate()` VPC parity).
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
