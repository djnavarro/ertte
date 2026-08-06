# ertte design history

This file is a condensed historical record of completed design
decisions: what was tried, what was rejected, and why. It exists for
context in future sessions, not as a changelog or PR log. Current-state
facts that came out of this history (what the API looks like today)
live in `AGENTS.md`, not here.

## API naming: `ertte_model()` split into `ertte_aft()`/`ertte_coxph()`

The design issue (#1) allowed for an optional `coxph()` semi-parametric
engine alongside the primary `survreg()`-based AFT engine. Rather than
one constructor with a `dist`/`engine` argument, they became two
separate, engine-specific constructors: Cox PH is structurally
different enough (semi-parametric, no location-scale structure, no
`dist` argument) that folding it into one constructor would be
misleading. There's no back-compat alias for the old `ertte_model()`
name -- the package was early enough in development that a clean
rename was preferred over a shim.

Downstream generics kept single shared names
(`ertte_predict()`/`ertte_fun()`), dispatching per engine via a shared
`"ertte_model"` superclass ahead of an engine-specific subclass
(`"ertte_aft"`/`"ertte_coxph"`) -- mirroring the base-R idiom of
`lm()`'s class `"lm"`. `simulate()` stayed a single
`simulate.ertte_model()` method for both engines, with per-engine
mechanics dispatched one level down via a new internal
`.ertte_simulate_draws()` S3 generic.

Later, `ertte_select_distribution()` was renamed to
`ertte_aft_select_distribution()` (2026-08-05) -- it was the one
AFT-exclusive function whose name didn't already say so (unlike
`ertte_aft()`/`ertte_coxph()` themselves, or the `dist` argument),
which read as though it might apply model-agnostically now that the
package has two engines. It never did, and never will, since Cox PH
has no `dist` to select over. No back-compat alias, matching the
`ertte_model()` precedent.

## Implementing the `ertte_coxph()` engine

Getting `survival::survfit()`/`basehaz()` to work at all on a
`coxph()` fit required fitting `ertte_coxph()` with `model = TRUE`:
`survfit.coxph()` otherwise tries to reconstruct the model frame by
re-evaluating `object$call$data`, which fails for the same reason
`stats::update()` fails on these fits (see "Refitting without
`stats::update()`" in `AGENTS.md`) -- the captured call refers to
`ertte_coxph()`'s own local `formula`/`data` bindings, not anything
visible in `survfit()`'s caller's frame. Storing the model frame
directly sidesteps that.

`ertte_predict.ertte_coxph()` uses `survival::survfit(object, newdata,
conf.int = conf_level)`: this computes each row's survival curve from
the fitted baseline hazard (Breslow/Efron, matching `object$method`)
and linear predictor, evaluates it at `time` via `summary(...,
extend = TRUE)` (letting `time` exceed the last observed follow-up,
holding survival constant beyond it), and takes confidence intervals
from `survfit()`'s own log-transform CI rather than a hand-rolled Wald
interval -- genuinely different from `ertte_predict.ertte_aft()`
methodologically, expected given the different model structure.

`ertte_fun.ertte_coxph()` evaluates `S(t | x) = S0(t)^exp((x -
xbar)'param)` via `survival::basehaz(object, centered = TRUE)` (held
constant beyond the last observed time) and `object$means` (the
covariate means `coxph()` centers on, needed because `basehaz()`'s
baseline is relative to that, not to `x = 0`). One genuine engine
difference surfaced here: Cox models have no intercept (it's absorbed
into the baseline hazard), but `stats::model.matrix()` on the model's
`terms()` adds one anyway since `terms()` doesn't record its absence
-- `ertte_fun.ertte_coxph()` drops that column explicitly so
`ncol(mm)` matches `length(coef(object))`.

`simulate()` support (`.ertte_simulate_draws.ertte_coxph()`) samples
event times by inverse-CDF sampling on the *cumulative hazard* scale
rather than a parametric quantile function: `S(t | x) = u` rearranges
to `H0(t) = -log(u) / exp(lp)`, so `.ertte_coxph_invert_basehaz()`
inverts the fitted (step-function) baseline cumulative hazard at that
target value. Verified by prototyping: simulating many draws at the
baseline covariate profile and comparing the empirical proportion
surviving past a given time to the fitted `S0(t)` reproduces it
closely.

## `ertte_add_term()`/`ertte_remove_term()` were AFT-hardcoded

They used to refit by calling `ertte_aft(...)` directly (not
`stats::update()`, which doesn't work here -- see above). Fixed by
introducing the internal `.ertte_refit()` S3 generic (dispatching on
`mod`'s class) so SCM works for `ertte_coxph` models too. Fixing this
surfaced a genuine engine difference: `stats::drop.terms()`'s output
(a `terms` object) can be passed straight back into `survreg()` as
`formula`, but `coxph()` rejects it (errors in
`terms.formula()`/`ExtractVars`) -- `ertte_remove_term()` now converts
it via `stats::formula()` first, which works for both engines.

## Power-function covariate parameterisation (`ertte_power()`)

The design issue called for "continuous covariates as power
functions, categorical covariates as factors," flagging CIs on the
transformed parameter (delta method vs profile likelihood) as an open
question. The companion `emaxnls` package was checked for precedent
first and turned out *not* to implement power functions either -- its
covariate model is plain additive/linear per structural parameter, fit
via `nls()` -- so this was new design, not a port.

The key insight: unlike `emaxnls`'s genuinely nonlinear structural
parameters, both `ertte_aft()` (log-location-scale AFT) and
`ertte_coxph()` (Cox PH) are already linear in their covariates on the
model's natural (log-time / log-hazard-ratio) scale. A power-function
effect -- `T = T_ref * (x/ref)^theta` (AFT) or `h(t|x) = h0(t) *
(x/ref)^theta` (Cox) -- is, after taking logs, exactly a linear term in
`log(x/ref)`. So `ertte_power(x, ref = NULL)` just returns
`log(x/ref)`: the fitted coefficient *is* the power exponent directly,
and its ordinary Wald CI already is the CI on that exponent -- no delta
method or profile likelihood needed, resolving the design issue's open
question for these two model families specifically.

`ref` defaults to `median(x, na.rm = TRUE)` (pop-PK/NONMEM convention);
every non-missing `x` must be strictly positive, which rules out using
it on covariates with a placebo/zero group (`dose`/`aucss`/`cmaxss`) --
by design, since the power-function language targets the *covariate
model* (age, weight, ...), not the primary exposure metric.
Predict-time consistency (reusing the original fit's `ref`) is handled
by a `makepredictcall.ertte_power()` method -- the same mechanism
`stats::poly()`/`splines::ns()` use -- which required no changes to any
model-matrix-building code elsewhere in the package, since all of it
already builds design matrices via `model.matrix()`/`predict()`/
`survfit()` against `object$terms`, and `model.frame()` already honours
a `terms` object's `predvars` attribute. `ertte_add_term()`/
`ertte_remove_term()`/SCM also needed no changes, since their term
handling already operates generically on formula term-labels. Nothing
prevents combining a plain linear term (`age`) and a power term
(`ertte_power(age)`) for the same variable -- deliberately left to the
user's judgement.

## Scalar E-R reductions: `ertte_landmark()` and `ertte_rmst()`

Phase 2 of the design issue (scalar E-R views of a TTE endpoint,
reusing erplots' `er_plot()`/`er_vpc()` grammars) shipped in stages.

**Landmark-binary** (`ertte_landmark()`) came first: `P(event by t*) =
1 - S(t*)` is a decreasing monotonic transform of `ertte_predict()`'s
own survival-probability output, so the confidence interval bounds
simply swap -- no recomputation, no new edge-case handling, since
whatever validity `ertte_predict()`'s interval has for a given engine
carries through unchanged. It's deliberately *not* a generic (unlike
`ertte_predict()`/`ertte_fun()`) -- it needs no engine-specific logic
of its own, since it delegates entirely to `ertte_predict()`, which
already dispatches on the subclass.

`er_predict.ertte_model()` forwards to it, threading a required
`landmark_time` argument through `...` (erplots' fixed `er_predict()`
signature has no TTE-specific argument slot).

**RMST** (`ertte_rmst()`) came next, as a genuine generic (unlike
`ertte_landmark()`), since computing an area under the curve needs the
whole survival curve, not a single `ertte_predict()` call.

- *AFT method*: `fit_rmst` is `stats::integrate()` of the closed-form
  `S(t|x)` from 0 to `tau`; `se_rmst` is an analytic delta method
  differentiating under the integral sign, propagating only `Var(mu)`
  -- the same simplification `ertte_predict.ertte_aft()` already
  makes. `.ertte_dist_info()` gained a `dbase` (density) entry for
  this gradient.
- *Cox method*: since the fitted baseline hazard is a right-continuous
  step function, `fit_rmst` is an *exact* finite sum of rectangle
  areas up to `tau`, not a numerical-quadrature approximation. The
  genuinely hard part was `se_rmst`, not the integral.

**Deriving the Cox `se_rmst` -- a real methodological finding.** The
first candidate was reusing `survival:::survmean()` (behind
`print.survfit(fit, rmean = tau)`) on a `survfit(coxph_object, newdata
= ...)` object. This works mechanically and was reimplemented against
public `survfit()` fields only (`time`, `surv`, `n.risk`, `n.event`) to
avoid the unexported function -- but reading `survmean()`'s source
revealed it estimates the variance-of-mean via a Greenwood-type term
computed from `n.risk`/`n.event`, which are **shared identically
across every covariate profile** (confirmed via `identical()`) -- i.e.
population-level risk sets from the shared baseline hazard, not
profile-specific. It never touches `survfit()`'s own `std.err` field,
which *does* correctly combine coefficient and baseline-hazard
uncertainty per profile (confirmed: `sf$logse == TRUE` and
`sf$cumhaz == -log(sf$surv)` exactly, so `std.err(t)` genuinely
estimates `SE[H(t|x)]`). For an extreme/leveraged covariate profile,
`survmean()`'s SE was found to be roughly 15x *smaller* than a
300-replicate nonparametric bootstrap SE for the same quantity -- a
real, not cosmetic, understatement of uncertainty.

Fixed by keeping `survmean()`'s exact rectangle/tail-weighted-sum
construction (`.ertte_rmst_pfun_delta()`) but substituting the
variance-increment source: `diff(std.err^2)` (the increment of the
profile-specific `Var[H(t|x)]`) in place of the population Greenwood
term. This is an approximation -- the coefficient-uncertainty
component of `H(t|x)` is really a single random direction shared
across every `t`, not a sum of independent per-jump increments -- but
cross-validated against the same bootstrap across two contrasting
covariate profiles, it tracked the bootstrap SE substantially more
closely than both `survmean()`'s naive population term and a cheaper
alternative that holds the baseline hazard fixed and only redraws
coefficients from `MVN(coef(object), vcov(object))` (which
over/understates uncertainty inconsistently depending on how extreme
the covariate profile is). See `PLAN.md` for the residual caveat.

CI construction on both engines is a symmetric Wald interval on the
RMST scale, *not* automatically bounded to `[0, tau]` the way
`ertte_predict()`'s survival-probability intervals are bounded to
`[0, 1]` by construction -- a known limitation, see `PLAN.md`.

**Extreme-`tau` numerical instability (issue #12).** Investigating a
suspected raw `integrate()` error for absurd `tau` found something
worse: for sufficiently large `tau`, `stats::integrate(s_fun, 0, tau,
mu = mu)`'s adaptive quadrature silently returned `0` -- no error, no
warning -- once `[0, tau]` was many orders of magnitude wider than the
region where the survival curve is non-negligible (confirmed: broke
between `tau = 1e5` and `5e5` for a representative Weibull fit). This
violated RMST's own monotonicity with no signal to the caller. Fixed
by reparameterising both integrals in `ertte_rmst.ertte_aft()` onto the
`u = log(t)` scale (`t = exp(u)`, `dt = exp(u) du`): since the upper
bound becomes `log(tau)`, the integration domain only grows
logarithmically with `tau`, keeping quadrature reliable up to at least
`tau ~ 1e9` for the same model/profile -- verified across all four
supported distributions, and cross-checked against a known
mathematical identity (the gradient integral converges to `fit_rmst`
itself as `tau -> infinity`). Even after this fix, `tau` several more
orders of magnitude out was found to occasionally destabilise the
delta-method gradient specifically (a much narrower, "bump"-shaped
integrand than `fit_rmst`'s broad-plateau curve, harder for adaptive
quadrature once the domain is stretched far enough) -- rather than
chase a fully bulletproof numerical guarantee for a pathological input
space, `ertte_rmst.ertte_aft()` now warns
(`.ertte_check_extreme_aft_tau()`) if any `tau` exceeds 10,000x the
fitting data's last observed follow-up time. `ertte_predict.ertte_aft()`/
`ertte_fun.ertte_aft()` were checked and found *not* to share this
failure mode, since neither ever calls `stats::integrate()` -- they
evaluate the base distribution's CDF/density directly, which saturates
cleanly to `0`/`1`.

**VPC scalar-reduction parity.** `er_simulate.ertte_model()`'s old
per-row `sim_time`/`sim_event` shape didn't satisfy the landmark/RMST
VPC's needs (the design issue notes a landmark-binary VPC "likewise
reuses `er_vpc()` unchanged"). Fixed via a new internal
`.ertte_simulate_scalar_resp()`: `sim_resp` is built from each
replicate's already-censored `sim_time`/`sim_event` (a genuinely
ambiguous outcome relative to `landmark_time`/`tau` becomes `NA`,
handled correctly by `er_vpc_add_simulated()`'s `mean(..., na.rm =
TRUE)`); `fit_resp` reuses `ertte_fun(object)` evaluated at each
replicate's own sampled coefficient draw (recovered from the `coef_*`
columns `.ertte_simulate_draws()` already attaches). For RMST, since
there's no single evaluation point, `.ertte_rmst_fit_resp_curve()`
evaluates `ertte_fun()` on a fixed 64-point grid and applies the
composite trapezoidal rule -- a deliberately coarser approximation than
`ertte_rmst()`'s own point estimate, acceptable since `fit_resp` here
is illustrative, not something a CI is built from. Verified
empirically: for a single covariate profile with `nsim = 3000`,
`mean(sim_resp, na.rm = TRUE)` and `mean(fit_resp)` both landed within
roughly 1-2% of `ertte_landmark()`'s/`ertte_rmst()`'s own fitted value,
for both engines. The IPCW/pseudo-value alternative to the
complete-case convention was investigated and rejected -- see
`PLAN.md`.

## erplots#10/#11: `predict_args`/`summary_args`/`simulate_args`

End-to-end testing of the scalar E-R reductions against erplots found
`er_plot_add_model(mod, landmark_time = 90)` erroring even though
`er_predict.ertte_model()`'s contract worked correctly when called
directly: `er_plot_add_model()`'s `...` was captured only for its
*style builder*, never forwarded to `er_predict()` itself. The
identical gap affected `er_plot_add_summary()`/`er_summary()` (which
forwarded *no* arguments at all, not even `conf_level`) and
`er_vpc_add_simulated()`/`er_simulate()`. Filed as
[erplots#10](https://github.com/djnavarro/erplots/issues/10), fixed
upstream by [erplots#11](https://github.com/djnavarro/erplots/pull/11)
(merged 2026-08-04): dedicated `predict_args`/`summary_args`/
`simulate_args = list()` arguments, spliced into the corresponding
generic call via `rlang::exec()`, keeping `...` exclusively for the
style builder. This needed **no ertte-side code changes** -- ertte's
methods already declared `...`, and once erplots started splicing
those named lists in, the existing methods picked it up automatically.

## `er_predict_survival.ertte_model()` (issue #13)

erplots' `er_tte_add_model()` (a parametric `S(t)` curve/ribbon overlay
on a Kaplan-Meier plot) needed a fifth interoperability generic beyond
`er_predict()`/`er_simulate()`/`er_summary()`:
`er_predict_survival(model, newdata, time_grid, conf_level, ...)`.
Unlike `er_predict()`, `newdata` here carries only covariate profiles
(no time column); `time_grid` is a separate numeric vector crossed
against `newdata` inside the method -- the contract's return shape was
deliberately designed by erplots to mirror `ertte_predict()`'s own.

`er_predict_survival.ertte_model()` is consequently close to a direct
pass-through to `ertte_predict()`, which already returns exactly that
column set/row order for both engines. The one wrinkle:
`er_tte_add_model()`'s default `time_grid` spans `object$time$limits`,
whose lower end is `0` (the conventional KM origin, `S(0) = 1`) -- but
`ertte_predict()` rejects a non-positive `time` outright (`log(time)`/
the baseline-hazard lookup are undefined there for either engine).
Rather than loosen `ertte_predict()`'s own contract for this one
caller, `er_predict_survival.ertte_model()` special-cases `time_grid`
entries of exactly `0` (returning `S(0) = 1` directly, no model
evaluation needed), calling `ertte_predict()` unmodified for every
strictly positive grid point and interleaving the two back into
`time_grid`'s original per-profile order. A zero-row `newdata` is also
handled directly rather than relying on `ertte_predict()`'s own
zero-row handling, since the interleaving logic needs `nrow(newdata)`
up front regardless.

Verified end-to-end against `erplots::er_tte()`/`er_tte_add_curve()`/
`er_tte_add_model()`, including the stratified case and the
quantile-binned-numeric-strata case. `ertte_coxph()`'s existing
all-censored/single-level-factor edge cases (see "Stress-test
findings" below) were checked against this method and no new gaps were
found -- since the method is a thin pass-through for every strictly
positive `time_grid` entry, existing per-engine guards carry through
unchanged. One genuine asymmetry is documented rather than silently
relied upon: a `time_grid` of exactly `0` bypasses every guard,
returning `S(0) = 1` even for an otherwise-invalid model -- correct
(a trivially true fact independent of model validity), but means an
all-censored `ertte_coxph` model can successfully answer `time_grid =
0` even though every other value errors. In practice this rarely
matters, since `er_tte_add_model()`'s default `time_grid` always spans
`0` and positive values, so the guard still fires for the realistic
default case.

## Administrative-censoring simulation refinement

`.ertte_simulate_draws()` used to cap every simulated event time at
that row's *observed* exit time, whether that row was itself an event
or a censoring. Splitting by the observed `event` indicator sharpened
the diagnosis: for censored rows, the observed exit time genuinely
*is* the true censoring time; for event rows, it isn't -- it's when the
event happened, not the (necessarily later, unobserved) administrative
censoring horizon -- so capping there leaked the observed outcome into
the simulation and biased simulated-vs-observed comparisons (e.g. a
VPC) toward looking more similar than the fitted model actually
implies.

`simulate.ertte_model()`/`.ertte_simulate_draws()` gained an optional
`censor_time` argument (a genuine per-row administrative follow-up
time, validated by `.ertte_check_censor_time()`); when supplied, it
caps every row uniformly regardless of observed event status. Absent
it, a new shared `.ertte_apply_admin_censoring()` helper caps censored
rows at their observed exit time (unchanged) but leaves event rows
uncensored -- a genuine, confirmed-with-the-maintainer default-behaviour
change, still an approximation for event rows, just a less biased one.
`ertte_data` gained a demonstration column, `admin_censor` (fixed at
180 for every row, reflecting the study-wide cutoff `.make_ertte_data()`
already used internally but didn't retain) -- confirmed this didn't
perturb any other column's simulated values under the same `seed =
111L`.

One real bug surfaced implementing the Cox engine's version:
`.ertte_coxph_invert_basehaz()` returns `Inf` whenever a simulated
draw would need to survive past the fitted baseline hazard's support
to "fail" -- previously always silently absorbed by the old blanket
`pmin(sim_time_raw, obs_time)` cap, but left as a bare `Inf`
(`sim_event = 1`, an event at time infinity) once event rows stopped
being capped by default. Fixed by substituting `max(bh$time)` for
`Inf` before applying `censor_time`/`obs_time`, then forcing
`sim_event = 0` for any row that hit this extrapolation boundary --
consistent with the flat baseline-hazard extrapolation
`ertte_predict.ertte_coxph()` already uses. The AFT engine's equivalent
failure mode (`qbase()` at `u` near 1) was left unguarded -- the risk
is negligible there (a continuous quantile function, vs. the Cox
engine's empirical/step-function baseline hazard, which realistically
exhausts its support whenever `nsim` is reasonably large).

## SCM selection criteria: p-value, AIC, BIC

`ertte_scm_forward()`/`ertte_scm_backward()` gained a `criterion`
argument (`"p-value"` (default), `"aic"`, or `"bic"`), mirroring the
companion `emaxnls` package's own development version rather than
introducing a new design. `threshold` is used only for `"p-value"`
mode and silently ignored for `"aic"`/`"bic"`. In IC mode, the internal
per-step helpers compare each candidate's refit against a running
`best_metric` (the current model's AIC/BIC, updated as better
candidates are found), adding/removing a term only if doing so
strictly decreases the chosen IC. Since a likelihood-ratio p-value
plays no role in IC-based selection, it's simply not computed under
`"aic"`/`"bic"` -- `term_p_value` is left `NA` in the history for that
step, matching `emaxnls`'s own behaviour. `model_aic`/`model_bic` are
always recorded regardless of which criterion drove selection.
`ertte_scm_history()` gained a `criterion` column (`NA` for the
base-model row and any pre-existing history predating this feature).

## Stress-test findings and fixes

Two rounds of adversarial edge-case testing surfaced ten genuine
issues (#3-#12), all now fixed and closed.

**Round 1 (#3-#7):**

- **#7 -- `ertte_fun()` performed no `time` validation**, unlike
  `ertte_predict()`. A new shared `.ertte_check_time()` helper
  (consolidating logic previously duplicated inline in both
  `ertte_predict()` methods) is now also called at the top of both
  `ertte_fun()` methods' returned closures. Before this fix, a
  non-positive `time` returned silently: `NaN` for AFT (via
  `log(time)`), or `1` for Cox (as if survival were guaranteed).
- **#6 -- `ertte_fun()`'s closure argument order (`param` before
  `data`)** was inconsistent with the rest of the API. Both closures
  are now `function(data = NULL, time, param = NULL)`, matching the
  "data-first" convention `ertte_predict()`/`ertte_landmark()`/
  `ertte_rmst()` already use. No caller inside the package used
  positional arguments, so this was a pure signature change with no
  follow-on breakage.
- **#5 -- `ertte_coxph()` didn't validate `time > 0`.** A new
  `.ertte_check_response_time()` helper is called at the top of
  `ertte_coxph()`, before `survival::coxph()` is invoked.
  `ertte_aft()` was left unchanged: it already gets equivalent
  validation "for free" via `survreg()`'s own internal check
  (`log(time)` is undefined for non-positive values) -- there was no
  gap to close there, just an engine inconsistency this fix resolves
  from the Cox side.
- **#4 -- a cryptic single-row error from `ertte_coxph()`.**
  `survival::coxph()` crashes with `'x' must be an array of at least
  two dimensions` (from `rowSums()` on its own post-fit influence
  diagnostics, which degrades from a matrix to a plain vector with
  only one observation) when fitting on exactly one usable row. A new
  `.ertte_check_coxph_data_size()` helper checks up front (via
  `stats::model.frame()`, not raw `nrow(data)`, since a raw-row count
  can be insufficient once missing covariates are accounted for) that
  at least 2 usable rows are available.
- **#3 -- empty `candidates` in `ertte_aft_select_distribution()`**
  silently fit no candidates and returned a degenerate result instead
  of erroring. A new `.ertte_check_dist_candidates()` helper validates
  up front.

A follow-up pass found two consolidation opportunities (no behaviour
change): `.ertte_check_dist_candidates()` and the pre-existing
`.ertte_check_candidates()` both started with an identical "non-empty
character vector, no missing values" check, factored out into a shared
`.ertte_check_nonempty_character()`. Separately, the new
`.ertte_surv_vars_from_formula()` (#5) duplicated the pre-existing
`.ertte_response_vars()`'s `Surv()`-parsing logic; `.ertte_response_vars()`
now delegates to it.

**Round 2 (#8-#12):**

- **#8/#9 -- SCM misattributed or crashed on a bad candidate.**
  `.ertte_once_forward()`/`.ertte_once_backward()` now wrap each
  candidate's refit in `tryCatch()`. A candidate whose refit throws
  (e.g. a single-level factor triggering a `contrasts` error) is now
  skipped with a warning quoting the actual error, rather than
  crashing the entire search (#9). A candidate that can't be
  added/removed at all (e.g. references a variable not in the fitting
  data) is now detected directly by comparing the refit's formula to
  the current model's, rather than misdiagnosing it via an `anova()`
  comparison against itself that produced a misreported "aliased/
  collinear" `NA` p-value (#8).
- **#10 -- zero-row `newdata` handling was inconsistent across
  engines.** `ertte_aft()`'s methods already worked by incidental
  behaviour; `ertte_coxph()`'s didn't (`survival::survfit()` errors on
  a zero-row `newdata` with a cryptic message). Both engines now
  short-circuit to a zero-row tibble with the expected columns before
  calling into `survival`.
- **#11 -- `conf_level = 0`/`1` worked for AFT but not Cox.**
  `ertte_aft()`'s methods already handle these degenerate endpoints
  correctly via their CDF back-transform, but `survival::survfit()`'s
  own `conf.int` argument rejects exactly 0 or 1. Fixed per-function,
  matched to how each uses `survfit()`'s output:
  `ertte_predict.ertte_coxph()` special-cases 0/1 directly;
  `ertte_rmst.ertte_coxph()` (which never uses `survfit()`'s own
  `$lower`/`$upper`) always passes a fixed valid placeholder and lets
  its existing `z_scale` logic handle 0/1 naturally.
- **#12** -- the `ertte_rmst.ertte_aft()` extreme-`tau` numerical
  instability; see "Scalar E-R reductions" above.

**Edge-case coverage (all-censored data, single stratum, heavy ties)**
was also added as regression tests, surfacing:

- **All-censored data (0 events) broke `ertte_coxph()`'s downstream
  methods.** `ertte_aft()` degrades gracefully (a `survreg()` "did not
  converge" warning, but still returns sensible predictions).
  `ertte_coxph()` fit fine, but `ertte_predict()`/`ertte_fun()`/
  `simulate()` used to fail with a cryptic `model.frame.default()`
  error -- `coxph(model = TRUE)` still leaves `object$model` unset with
  zero events, defeating the `model = TRUE` workaround. Fixed with a
  new `.ertte_check_coxph_nevent()` guard, called at the top of the
  three baseline-hazard-based methods, aborting informatively instead.
  The constructor itself is left alone -- fitting on all-censored data
  is still allowed.
- **Single-level factor covariates** fit fine on both engines
  (aliased coefficient reported as `NA`), but downstream behaviour
  genuinely differs by *method*: `predict.survreg()` (AFT) propagates
  the `NA` through; `survival::survfit()` (Cox) silently drops the
  aliased column and returns finite predictions anyway;
  `ertte_fun.ertte_coxph()`'s manual `mm %*% param` propagates `NA`
  like the AFT method. None of this crashes and isn't treated as a bug
  -- pinned down with regression tests rather than left as a surprise.
- **Heavy ties in event times** caused no issues on either engine.

## Documentation restructuring

`AGENTS.md` originally accreted "used to be X, now Y" narrative and
scoped-out future work directly, growing past 900 lines. Split
following the pattern established in the companion `erplots` repo:
`AGENTS.md` kept as a lean, current-state architecture reference;
completed design decisions and their rationale moved to
`.agents/HISTORY.md` (this file); scoped-out future work moved to
`.agents/PLAN.md`. Both new files are excluded from the built package
via `.Rbuildignore` (`^\.agents$`), matching how `AGENTS.md` itself is
already excluded; pkgdown only scans the package root (and `.github/`)
for stray `*.md` pages, so no change to `tools/pkgdown-postbuild.R` was
needed for the new subdirectory.
