# AGENTS.md

## What this package is

ertte is the time-to-event (TTE) member of the exposure-response (E-R)
package family, alongside [erglm](https://github.com/djnavarro/erglm)
(GLM E-R models) and [emaxnls](https://github.com/djnavarro/emaxnls)
(Emax / logistic-Emax models), all of which plug into
[erplots](https://github.com/djnavarro/erplots) for visualisation. See
[GitHub issue #1](https://github.com/djnavarro/ertte/issues/1) for the
full design/scoping discussion this package implements.

This first pass covered the issue's suggested phasing step 1 only:
**core modelling, methods, SCM, and simulation -- no plotting code**.
Two later phases from the issue were originally deferred; Phase 2 is now
fully implemented, Phase 3 is not started (see "Planned work" below):

- Phase 2: scalar E-R views of a TTE endpoint (landmark-binary / RMST)
  reusing erplots' existing `er_plot()`/`er_vpc()` grammars -- **done**,
  including both scalar reductions' point/CI predictions
  (`er_predict.ertte_model()`) and VPC simulation
  (`er_simulate.ertte_model()`).
- Phase 3: a new `er_tte()` KM/survival-curve plotting grammar in the
  separate `erplots` repo -- not started.

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
- **`er_predict.ertte_model()`/`er_simulate.ertte_model()`'s real
  contract -- landmark-binary, RMST, and their VPC parity are all now
  implemented.** Phase 2's Workstream B1 (scalar E-R views of a TTE
  endpoint) is now fully done:
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

  `er_predict.ertte_model()` now also forwards to `ertte_rmst()` when a
  `tau` argument (rather than `landmark_time`) arrives through `...`,
  since RMST is the other scalar E-R reduction the design issue
  mentions. Supplying both `landmark_time` and `tau` in the same call
  errors (they select different reductions), and -- unlike
  `ertte_rmst()`'s own vectorised `tau` argument, which can evaluate
  several horizons at once -- only a single `tau` value is accepted
  here, since erplots' scalar E-R grammar expects exactly one row per
  `newdata` row. Both branches rename their reduction's fitted-value
  column (`ertte_landmark()`'s `fit_resp` already matches; `ertte_rmst()`'s
  `fit_rmst` is renamed to `fit_resp`) to the single shared name
  erplots' plotting grammar expects, regardless of which reduction
  produced it.

  **`er_simulate.ertte_model()` scalar-VPC parity -- now implemented,
  for both landmark-binary and RMST.** The design issue notes a
  landmark-binary VPC "likewise reuses `er_vpc()` unchanged", which
  needs `er_simulate.ertte_model()` to return landmark-transformed
  draws (`fit_resp`/`sim_resp`), not its previous per-row
  `sim_time`/`sim_event` shape. Supplying `landmark_time` or `tau`
  through `...` (mutually exclusive, matching `er_predict.ertte_model()`'s
  scheme above) now does this, via a new internal
  `.ertte_simulate_scalar_resp()` (in `R/er-methods.R`):
  - **`sim_resp`** (erplots' hard requirement for `er_vpc_add_simulated()`,
    per `?erplots::er_model_interface`) is built directly from each
    replicate's *already-censored* `sim_time`/`sim_event` -- reusing the
    same administrative-censoring convention `.ertte_simulate_draws()`
    applies elsewhere (`.ertte_apply_admin_censoring()`), which is
    exactly what makes the simulated data comparable to a real,
    similarly-censored observed study, the entire point of a VPC. A
    replicate whose simulated outcome is genuinely ambiguous relative to
    `landmark_time`/`tau` (censored strictly before it) becomes `NA` --
    the same complete-case convention a landmark/RMST analysis already
    has to apply to genuinely censored *observed* data (e.g. the manual
    `case_when()` construction `test-er-methods.R`'s landmark test uses
    to build its observed-side data), and one `er_vpc_add_simulated()`'s
    `mean(..., na.rm = TRUE)` aggregation handles correctly by simply
    excluding it. For RMST, the per-replicate individual quantity is
    `min(sim_time, tau)` when the outcome relative to `tau` is known (an
    event, or survival to/past `tau`) -- the same construction that
    gives `E[min(T, tau)] = RMST(tau)` in the population-level
    formalism (see `vignettes/articles/rmst.Rmd`).
  - **`fit_resp`** (optional, for erplots' spaghetti-style plots) reuses
    `ertte_fun(object)` (already implemented, engine-agnostic) evaluated
    at each replicate's own sampled coefficient draw -- recovered from
    the `coef_*` columns `.ertte_simulate_draws()` already attaches per
    replicate, so no new coefficient sampling was needed. For
    `landmark_time` this is a single `ertte_fun()` call per replicate
    (vectorised across every `newdata` row in that replicate at once).
    For `tau` (RMST), there's no single evaluation point -- the whole
    curve from 0 to `tau` needs integrating -- so a new
    `.ertte_rmst_fit_resp_curve()` evaluates `ertte_fun()` on a fixed
    64-point grid and applies the composite trapezoidal rule, still
    vectorised across every `newdata` row in the replicate at once.
    This is a deliberately coarser approximation than `ertte_rmst()`'s
    own point estimate (an exact step-function sum for the Cox engine,
    adaptive quadrature for AFT) -- an acceptable trade-off since
    `fit_resp` here is an illustrative, optional quantity, not something
    a confidence interval is built from.
  - Verified both point estimates empirically: for a single covariate
    profile with a large `nsim` (3000), `mean(sim_resp, na.rm = TRUE)`
    and `mean(fit_resp)` both landed within roughly 1-2% of
    `ertte_landmark()`'s/`ertte_rmst()`'s own fitted value for the same
    profile, for both engines.
  - Verified end-to-end against erplots' actual `er_vpc()`/
    `er_vpc_add_observed()`/`er_vpc_add_simulated(simulate_args =
    list(landmark_time = ...))` (and the `tau` equivalent) --
    `er_vpc_build()` renders a `ggplot` object with no ertte-side
    surprises, reusing the `predict_args`/`simulate_args` splicing
    erplots#11 already added (see below) with no further erplots-side
    changes needed.
  - A shared `.ertte_check_single_tau()` (in `R/utils-helpers.R`)
    validates `tau` for both `er_predict.ertte_model()` and
    `er_simulate.ertte_model()`, replacing the inline check
    `er_predict.ertte_model()` previously had (no behaviour change,
    just avoiding duplicating the same validation message twice).
  - **The complete-case (NA-for-ambiguous-censoring) convention was
    considered for a configurable IPCW alternative, and rejected for
    now.** Inverse-probability-of-censoring weighting (Robins-style:
    keep every known-outcome replicate but upweight it by `1/Ĝ(t)`,
    `Ĝ` a KM fit on the censoring indicator, to compensate for peers
    dropped as ambiguous) is the standard fix for informative censoring
    bias -- but checking erplots' actual aggregation code
    (`R/er-vpc-layer.R` in the `erplots` repo) confirmed every
    observed/simulated summary is a plain unweighted `mean(...,
    na.rm = TRUE)`; there's no weight column anywhere in the VPC
    contract for a per-replicate weight to reach. A pseudo-observations
    alternative (Andersen-Perme-style: reduce each censored outcome to
    a single continuous value compatible with ordinary averaging,
    avoiding the missing-weight-column problem) would fit the contract
    mechanically, but needs a leave-one-out KM/RMST jackknife
    recomputation across the sample -- doing that *per simulated
    replicate* (`nsim` separate leave-one-out passes over `n`
    observations) is a real computational cost, and conceptually
    murky besides: pseudo-values correct for censoring bias in an
    unknown population, but every simulated replicate here is already a
    fully known draw from the fitted model, so it's unclear what
    "correcting for censoring bias" would even mean applied to data
    ertte generated itself. Net conclusion: a real IPCW/pseudo-value
    toggle would need an erplots-side contract change first (a weighted
    aggregation path, mirroring the erplots#10/#11 pattern below) before
    it could be usefully added here. In the meantime, `censor_time` (via
    `simulate_args`) is the existing lever for reducing how much
    ambiguity arises in the first place -- supplying the true
    administrative follow-up horizon per subject, rather than falling
    back to the default (censored rows capped at their observed exit
    time, event rows left uncensored).

  **Caveat discovered end-to-end testing this against erplots'
  `er_plot()` -- now resolved upstream.** `er_plot_add_model(mod,
  landmark_time = 90)` used to error, even though
  `er_predict.ertte_model()`'s contract (above) was implemented
  correctly and worked when called directly: `er_plot_add_model()`'s
  `...` was captured only for its *style builder* (`config$dots`, per
  `?er_style`), never forwarded to `er_predict()` itself. The identical
  gap affected `er_plot_add_summary()`/`er_summary()` (which forwarded
  *no* arguments at all, not even `conf_level`) and
  `er_vpc_add_simulated()`/`er_simulate()`. Filed as
  [erplots#10](https://github.com/djnavarro/erplots/issues/10), this
  was fixed upstream by [erplots#11](https://github.com/djnavarro/erplots/pull/11)
  (merged 2026-08-04), which took the "Design B" option sketched on the
  issue: dedicated `predict_args`/`summary_args`/`simulate_args = list()`
  arguments on `er_plot_add_model()`/`er_plot_add_summary()`/
  `er_vpc_add_simulated()` respectively, spliced into the corresponding
  generic call via `rlang::exec(er_predict, ..., !!!predict_args)` (and
  analogously for the other two) -- keeping `...` exclusively for the
  style builder, with no ambiguity about which consumer a given named
  argument reaches.

  Confirmed this needed **no ertte-side code changes at all**: ertte's
  `er_predict.ertte_model()`/`er_summary.ertte_model()`/
  `er_simulate.ertte_model()` already declared `...` (to forward
  `landmark_time`/`censor_time` etc.), and erplots' `er_predict()`/
  `er_summary()`/`er_simulate()` generics already declared `...` too --
  once erplots started actually splicing `predict_args`/`summary_args`/
  `simulate_args` into those generic calls, the existing ertte methods
  picked it up automatically. Verified directly:
  `er_plot_add_model(mod, predict_args = list(landmark_time = 90))`,
  `er_plot_add_summary(model = mod, conf_level = 0.9)` (now genuinely
  reaches `er_summary.ertte_model()`'s `conf_level`, previously
  impossible), and `er_simulate(mod, newdata, simulate_args =
  list(censor_time = 200))`-style calls (via `simulate_args` on
  `er_vpc_add_simulated()`) all now work as intended, with the
  `er_plot()` pipeline rendering correctly end-to-end. A regression
  test (`test-er-methods.R`) exercises the `er_plot_add_model(mod,
  predict_args = list(landmark_time = 90))` path directly against
  erplots, skipping gracefully if the installed erplots predates
  `predict_args`. The `new_landmark_model()` wrapper-object workaround
  documented in earlier revisions of this file is no longer needed and
  has been removed from here; use `predict_args`/`summary_args`/
  `simulate_args` instead.
- **RMST (restricted mean survival time) -- now implemented**, via a
  new exported `ertte_rmst(object, newdata, tau, conf_level)` (in
  `R/ertte-rmst.R`). Unlike `ertte_landmark()`, this is a genuine
  generic (`ertte_rmst.ertte_aft()`/`ertte_rmst.ertte_coxph()`), since
  computing an area under the curve needs the whole survival curve, not
  a single `ertte_predict()` call at one time point -- it doesn't
  delegate to `ertte_predict()` the way `ertte_landmark()` does.
  - **AFT method**: `fit_rmst` is `stats::integrate()` of the closed-form
    `S(t|x)` from 0 to `tau`; `se_rmst` is an analytic delta method
    differentiating under the integral sign (`d/dmu RMST = integral of
    dbase(z)/scale`), propagating only `Var(mu)` from `predict(object,
    newdata, type = "linear", se.fit = TRUE)` -- the same
    ignore-`scale`-uncertainty simplification `ertte_predict.ertte_aft()`
    already makes. `.ertte_dist_info()` (in `R/utils-helpers.R`) gained a
    `dbase` (base-distribution density) entry alongside its existing
    `pbase`/`qbase`, needed for this gradient.
  - **Cox method**: since the fitted baseline hazard (and every
    covariate-adjusted curve) is a right-continuous step function,
    `fit_rmst` is an *exact* finite sum of rectangle areas between jump
    times up to `tau` -- not a numerical-quadrature approximation, contra
    this file's earlier framing of RMST as needing "numerical
    integration" for the Cox engine. The genuinely hard part turned out
    to be `se_rmst`, not the integral itself.
  - **Deriving `se_rmst` for the Cox method -- a real methodological
    finding.** The first candidate investigated was reusing
    `survival:::survmean()` (the internal function behind
    `print.survfit(fit, rmean = tau)`) on a `survfit(coxph_object,
    newdata = ...)` object. This works mechanically (produces
    profile-varying `rmean`/`se(rmean)`, vectorises across multi-row
    `newdata`), and was reimplemented against public `survfit()` fields
    only (`time`, `surv`, `n.risk`, `n.event`) to avoid depending on the
    unexported `:::` function -- but reading `survmean()`'s source
    revealed it estimates the variance-of-mean via a Greenwood-type term
    (`n.event / (n.risk * (n.risk - n.event))`) computed from `n.risk`/
    `n.event`, which are **shared identically across every covariate
    profile** (confirmed via `identical()`) -- i.e. population-level risk
    sets from the shared baseline hazard, not profile-specific. It never
    touches `survfit()`'s own `std.err` field, which *does* correctly
    combine coefficient and baseline-hazard uncertainty per profile
    (confirmed empirically: `sf$logse == TRUE` and `sf$cumhaz ==
    -log(sf$surv)` exactly, so `std.err(t)` genuinely estimates
    `SE[H(t|x)]`). Net effect: for an extreme/leveraged covariate
    profile, `survmean()`'s SE was found to be roughly 15x *smaller*
    than a 300-replicate nonparametric bootstrap SE for the same
    quantity -- a real, not cosmetic, understatement of uncertainty.
    Fixed by keeping `survmean()`'s exact rectangle/tail-weighted-sum
    construction (now `.ertte_rmst_pfun_delta()` in `R/ertte-rmst.R`)
    but substituting the variance-increment source: `diff(std.err^2)`
    (the increment of the profile-specific `Var[H(t|x)]`) in place of
    the population Greenwood term. This is an approximation -- the
    coefficient-uncertainty component of `H(t|x)` is really a single
    random direction shared across every `t`, not a sum of independent
    per-jump increments, so treating its increments as accumulating
    independently (the same simplifying assumption the classic Greenwood
    formula makes) isn't exactly right -- but cross-validated against
    the same bootstrap across two contrasting covariate profiles, it
    tracked the bootstrap SE substantially more closely than both
    `survmean()`'s naive population term and a cheaper alternative that
    holds the baseline hazard fixed and only redraws coefficients from
    `MVN(coef(object), vcov(object))` (which over/understates uncertainty
    inconsistently depending on how extreme the covariate profile is).
  - **CI construction, both engines**: `ci_lower`/`ci_upper` are
    symmetric Wald intervals on the RMST scale (`fit_rmst +/- z *
    se_rmst`), *not* automatically bounded to `[0, tau]` the way
    `ertte_predict()`'s survival-probability intervals are bounded to
    `[0, 1]` by construction (their CDF back-transform keeps them there)
    -- a known, documented limitation, not yet addressed.
  - `ertte_rmst.ertte_coxph()` warns if any `tau` exceeds the last
    observed follow-up time across the fitted cohort: RMST integrates
    the *entire* curve up to `tau`, so the flat-baseline-hazard
    extrapolation convention `ertte_predict.ertte_coxph()` already uses
    for a single time point has a larger, more silent effect on an area.
    It also inherits `.ertte_check_coxph_nevent()`'s all-censored-data
    guard, same as `ertte_predict.ertte_coxph()`/`ertte_fun.ertte_coxph()`.
  - Tests live in `tests/testthat/test-ertte-rmst.R`, including a
    regression test pinning down that `fit_rmst` matches a hand-computed
    step-function integral, and a regression test confirming the derived
    Cox `se_rmst` is *not* numerically equal to `survival:::survmean()`'s
    naive population-level SE (i.e. that the fix above stays fixed).
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
    `.ertte_simulate_draws()` S3 generic (see `R/ertte-aft.R`).
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
  - `ertte_rmst()` -- a generic, with methods for both engines:
    `ertte_rmst.ertte_aft()` (quadrature + analytic delta method) and
    `ertte_rmst.ertte_coxph()` (exact step-function sum + a delta method
    on the profile-specific `std.err` field) -- see "Planned work" above
    for the derivation and the bootstrap cross-check behind the Cox
    method's confidence interval.

## Structure

- `R/ertte-aft.R` -- `ertte_aft()`, the `ertte_predict()`/`ertte_fun()`
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
- `R/ertte-rmst.R` -- `ertte_rmst()`: the restricted-mean-survival-time
  scalar E-R reduction, with `ertte_rmst.ertte_aft()` (quadrature +
  analytic delta method) and `ertte_rmst.ertte_coxph()` (exact
  step-function sum, via the internal `.ertte_rmst_pfun_delta()`
  helper, plus a delta method built on `survfit()`'s profile-specific
  `std.err` field) methods. Unlike `ertte_landmark()`, a genuine
  generic -- see "Planned work" above for why, and for the derivation
  of the Cox method's confidence interval. A user-facing explanation of
  the underlying formalism, assumptions, and caveats -- pitched at a
  pharmacometrician familiar with TTE models but not the statistical
  details of AFT/Cox PH -- lives in a website-only article,
  `vignettes/articles/rmst.Rmd` (see "Development workflow" below for
  why it's an article rather than a package vignette).
- `R/ertte-landmark.R` -- `ertte_landmark()`: the landmark-binary
  scalar E-R reduction (`P(event by t*) = 1 - S(t*)`) that
  `er_predict.ertte_model()` (`R/er-methods.R`) forwards to. A single
  function, not a generic -- it delegates entirely to
  `ertte_predict()`, which already dispatches per engine. See "Planned
  work" above for what's deferred (`er_simulate()` VPC parity). Has its
  own, shorter website-only article, `vignettes/articles/landmark.Rmd`
  (see "Development workflow" below), loosely modelled on
  `rmst.Rmd`'s structure.
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
  generics themselves. Has a website-only article,
  `vignettes/articles/scm.Rmd` (see "Development workflow" below),
  covering forward/backward selection, the LRT-based significance test,
  and how `ertte_power()` terms fit into a candidate set.
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

- **Longer explanatory write-ups (formalism, assumptions, caveats) go in
  website-only articles, not package vignettes.** `vignettes/articles/`
  (scaffolded via `usethis::use_article()`) holds plain `.Rmd` files with
  no `VignetteBuilder`/vignette-engine metadata -- `.Rbuildignore` excludes
  the whole directory from the built package, but pkgdown still discovers
  and renders them into an "Articles" section on the site. Used for
  `ertte_rmst()`'s methodology write-up (`vignettes/articles/rmst.Rmd`),
  since its confidence-interval derivation (see "Planned work" above) is
  too long and too caveat-heavy for `?ertte_rmst`'s `@details`, but doesn't
  need the R CMD build/check overhead (or `VignetteBuilder: knitr` in
  `DESCRIPTION`) a real package vignette would add. Render locally with
  `rmarkdown::render()` to check it knits before pushing -- pkgdown's build
  step doesn't run in this repo's CI, only on deploy.
  Three further articles now exist alongside `rmst.Rmd`, all in
  `vignettes/articles/` and registered in `_pkgdown.yml`'s `articles:`
  list (`overview`, `landmark`, `rmst`, `scm`, in that reading order):
  `overview.Rmd` (a `ertte_aft()`/`ertte_coxph()` primer, with a
  refresher on the survival/hazard functions and the AFT/proportional-
  hazards ideas each engine is built on -- the intended starting point
  for a new user, cross-linked from the other three), `landmark.Rmd`
  (`ertte_landmark()`, deliberately much shorter than `rmst.Rmd` since
  a single-time-point reduction has far less to explain -- delegates
  entirely to `ertte_predict()`, so most of its content is about the
  CI-swap trick and the landmark VPC's censoring convention rather than
  new numerical machinery), and `scm.Rmd` (`ertte_scm_forward()`/
  `ertte_scm_backward()`/`ertte_add_term()`/`ertte_remove_term()`,
  including how `ertte_power()` terms need no special handling in a
  candidate set). All three were verified to render cleanly via
  `rmarkdown::render()` against a locally reinstalled copy of the
  package (`devtools::install()`) -- worth remembering that rendering
  articles against a stale *installed* copy of ertte (rather than the
  current source) can silently exercise old behaviour or fail on
  recently-renamed functions, as happened once during this work.
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
