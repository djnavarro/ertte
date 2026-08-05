# Restricted mean survival time in ertte

``` r

library(ertte)
library(survival)
```

## Why RMST?

Most exposure-response summaries for a time-to-event (TTE) endpoint
focus on either a hazard ratio (from a Cox model) or a survival
probability at a fixed time (from
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)/[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)).
Both are useful, but neither directly answers a question that’s often
more clinically intuitive: *on average, how much event-free time does a
subject with a given exposure profile experience, up to some horizon?*

**Restricted mean survival time (RMST)** answers exactly that. For a
survival function $`S(t \mid x)`$ describing a subject with covariate
profile $`x`$, and a fixed horizon $`\tau`$, RMST is defined as

``` math
\text{RMST}(\tau \mid x) = \int_0^\tau S(t \mid x) \, dt.
```

Geometrically, this is the area under the survival curve between 0 and
$`\tau`$. It has units of time (e.g. “expected days event-free in the
first 90 days”), doesn’t rely on a proportional-hazards assumption to
interpret, and stays well-defined even when a study’s follow-up is too
short to observe the survival curve reach zero – which is exactly why
it’s “restricted” to $`\tau`$ rather than defined over $`[0, \infty)`$.

[`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
computes RMST (with a confidence interval) for a fitted
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
model, at one or more horizons `tau`, for one or more covariate profiles
in `newdata`. This article works through how it does that for each model
type, what assumptions the calculation relies on, and where the
confidence interval is more approximate than you might expect.

This assumes you already have some familiarity with fitting TTE models
via the `survival` package
([`survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)/[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html))
and with `ertte`’s two model engines; if you’re new to `ertte`, start
with
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)/[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)’s
documentation and
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
first.

## Two engines, two survival curves

`ertte` fits two kinds of TTE models, and
[`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
needs a genuinely different calculation for each, because they represent
$`S(t \mid x)`$ differently.

**[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)**
fits a parametric accelerated failure time (AFT) model via
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).
All four supported distributions share the same log-location-scale
structure,

``` math
\log(T) = \mu(x) + \sigma W,
```

where $`\mu(x)`$ is the linear predictor (intercept + covariate
effects), $`\sigma`$ is a scale parameter, and $`W`$ follows a fixed
“base” distribution (extreme-value for the Weibull/exponential case,
standard normal for log-normal, standard logistic for log-logistic).
This gives a fully parametric, closed-form survival function
$`S(t \mid x) = 1 - F\!\left(
\frac{\log t - \mu(x)}{\sigma}\right)`$, where $`F`$ is the base
distribution’s CDF.

**[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)**
fits a semi-parametric Cox proportional-hazards model via
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).
There’s no parametric form for the baseline hazard $`h_0(t)`$; it’s
estimated nonparametrically (Breslow’s estimator, the same one behind
[`survival::basehaz()`](https://rdrr.io/pkg/survival/man/basehaz.html)),
and the survival curve for a covariate profile $`x`$ is

``` math
S(t \mid x) = S_0(t)^{\exp\left((x - \bar x)'\beta\right)},
```

where $`S_0(t) = \exp(-H_0(t))`$ is the fitted baseline survival curve
and $`\bar x`$ is the mean covariate profile the partial likelihood was
centered on. Because $`H_0(t)`$ is estimated from a finite set of
observed event times, $`S(t \mid x)`$ is a **right-continuous step
function** – it only changes value at the observed event times in the
fitting data, however far $`x`$ is from $`\bar x`$.

That structural difference – closed-form vs. step-function – is what
drives most of the implementation differences below.

## Computing the point estimate

### AFT: numerical integration of a known curve

Since $`S(t \mid x)`$ has a closed form for the AFT engine,
[`ertte_rmst.ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
computes

``` math
\widehat{\text{RMST}}(\tau \mid x) = \int_0^\tau \hat S(t \mid x) \, dt
```

via [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html). This
is straightforward numerical quadrature on a smooth, fully known
function – there’s no approximation of consequence here beyond ordinary
quadrature error, which is negligible for a function this well-behaved.

### Cox: an exact sum, not an approximation

For the Cox engine, it’s tempting to assume the same kind of numerical
integration is needed, since the “curve” comes from an empirical
baseline hazard rather than a formula. In fact the opposite is true:
because $`\hat S(t \mid x)`$ is a **step function** with jumps only at
the observed event times $`t_{(1)} < t_{(2)} < \dots < t_{(m)}`$, the
area under it between two jumps is exactly a rectangle, and the whole
integral reduces to an *exact* finite sum:

``` math
\widehat{\text{RMST}}(\tau \mid x) = \sum_{k} \hat S\!\left(t_{(k)}^{-} \mid x\right) \cdot
\Big(\min(t_{(k+1)}, \tau) - t_{(k)}\Big).
```

[`ertte_rmst.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
computes this directly from `survival::survfit(object, newdata)`’s own
`time`/`surv` fields – there’s no
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) call
anywhere in this path, and no quadrature error to worry about. The only
genuine approximation in the point estimate shows up at the boundary,
discussed under “Extrapolating beyond the data” below.

## Confidence intervals: the harder half

Getting a point estimate is the easy part for both engines. Getting a
defensible confidence interval turned out to need real design decisions,
and the two engines’ methods are not analogous to each other.

### AFT: a delta method that ignores scale uncertainty

[`ertte_predict.ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
already reports a Wald confidence interval on the survival probability
by propagating the uncertainty in $`\mu(x)`$ (the linear predictor) and
*not* the uncertainty in $`\sigma`$ (the scale parameter) – a
simplification made throughout `ertte`, on the grounds that $`\sigma`$’s
sampling variability is typically small relative to $`\mu(x)`$’s and
that jointly propagating both would need a full delta-method gradient
rather than a simple `predict(..., se.fit = TRUE)` call.
[`ertte_rmst.ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
makes the identical simplification.

Differentiating the RMST integral with respect to $`\mu`$
(differentiating under the integral sign) gives

``` math
\frac{\partial}{\partial \mu} \text{RMST}(\tau \mid x) =
\int_0^\tau \frac{f\!\left(\frac{\log t - \mu(x)}{\sigma}\right)}{\sigma} \, dt,
```

where $`f`$ is the base distribution’s density.
[`ertte_rmst.ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
evaluates this gradient numerically (another
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) call),
then applies the standard delta-method variance,

``` math
\widehat{\text{Var}}\!\left[\widehat{\text{RMST}}(\tau \mid x)\right] \approx
\left(\frac{\partial \text{RMST}}{\partial \mu}\right)^{\!2} \widehat{\text{Var}}(\hat\mu),
```

with $`\widehat{\text{Var}}(\hat\mu)`$ coming from
`predict(object, newdata, type = "linear", se.fit = TRUE)`, exactly as
in
[`ertte_predict.ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_predict.md).
The resulting interval is a symmetric Wald interval, `fit_rmst`
$`\pm\ z \cdot`$`se_rmst`.

### Cox: why the obvious shortcut doesn’t work

For the Cox engine, the natural first idea is to reuse machinery
`survival` already ships: `print(survfit_object, rmean = tau)` (backed
by the unexported `survival:::survmean()`) already prints an
`rmean`/`se(rmean)` pair for a fitted `survfit` object. Calling this on
`survival::survfit(coxph_object, newdata = x)` does run, and does
produce a value that looks reasonable at first glance – it varies by
covariate profile, because the *survival curve itself* differs by
profile.

The problem surfaces on closer inspection. `survmean()`’s variance
calculation uses a Greenwood-type increment,

``` math
\hat h_k = \frac{d_k}{n_k (n_k - d_k)},
```

built entirely from `n.risk` ($`n_k`$) and `n.event` ($`d_k`$) – the
*number of subjects at risk and experiencing an event* at each observed
time, in the whole fitting cohort. These are properties of the shared
baseline hazard estimate, **not of the covariate profile being predicted
for**: two very different covariate profiles passed through the same
fitted Cox model get *identical* `n.risk`/`n.event` sequences, because
both curves are built from the same baseline hazard and the same
underlying risk sets. The formula never touches
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)‘s own
`std.err` field, which – for a `coxph`-based curve – genuinely does
combine both sources of uncertainty (the baseline hazard’s estimation
uncertainty *and* the regression coefficients’ estimation uncertainty)
correctly for that specific profile.

The practical consequence: for a covariate profile far from the average
(e.g. an unusually high exposure value), `survmean()`’s naive standard
error can be dramatically too small – in one check during development,
roughly 15-fold smaller than a 300-replicate nonparametric bootstrap
standard error for the same quantity. That’s not a subtle numerical
discrepancy; it’s a real, and potentially misleading, understatement of
uncertainty for exactly the covariate profiles an exposure-response
analysis often cares most about (e.g. a high-exposure arm).

### Cox: the corrected delta method

[`ertte_rmst.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
instead builds on
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)’s own
`std.err` field, which is on the cumulative-hazard scale: `std.err(t)`
is the estimated standard error of $`\hat H(t \mid x)`$, the
profile-specific cumulative hazard (confirmed by checking that
`sf$cumhaz` equals `-log(sf$surv)` exactly, and `sf$logse` is `TRUE`).
This *does* vary correctly by covariate profile, since it comes from the
same delta-method calculation `survfit.coxph()` uses internally to build
its own per-time confidence intervals.

The fix keeps `survmean()`’s overall construction – summing squared
“area remaining beyond $`t_k`$” terms, each weighted by a variance
increment at $`t_k`$ – but replaces the increment itself:

``` math
\hat h_k = \widehat{\text{Var}}\!\left[\hat H(t_{(k)} \mid x)\right] -
\widehat{\text{Var}}\!\left[\hat H(t_{(k-1)} \mid x)\right] =
\text{std.err}(t_{(k)})^2 - \text{std.err}(t_{(k-1)})^2.
```

This is a genuine improvement, but it’s still an **approximation**, and
it’s worth understanding exactly which assumption it makes. The classic
Greenwood-type construction (and this adapted version of it) implicitly
treats $`\hat H(t \mid x)`$ as accumulating via a sequence of
*statistically independent* increments over time – true, in an
appropriate asymptotic sense, for the martingale-based part of a
nonparametric hazard estimator. But for a Cox model, part of
$`\hat H(t \mid x)`$’s uncertainty comes from the regression
coefficients $`\hat\beta`$, and that part is really **one shared random
quantity affecting every time point together**, not something that
accumulates independently as $`t`$ increases. Treating its contribution
as if it did accumulate independently is a simplification, not an exact
result.

This was checked empirically during development: for two contrasting
covariate profiles, a 300-replicate nonparametric bootstrap (refitting
the Cox model on resampled data and recomputing RMST each time) was used
as a independent benchmark. The corrected delta-method standard error
tracked the bootstrap substantially more closely than either
`survmean()`’s naive version or a second, cheaper alternative that holds
the baseline hazard fixed and treats only $`\hat\beta`$ as random
(which, depending on how extreme the covariate profile is, can either
over- or under-state the true uncertainty). It was not, however, checked
across a wide range of sample sizes, censoring patterns, or covariate
profiles – treat it as a meaningfully better approximation, not a
proven-exact one.

## Extrapolating beyond the data

Both engines can, in principle, be asked for RMST at a horizon $`\tau`$
that exceeds the longest follow-up time observed in the fitting data.

- For
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md),
  the baseline hazard has no information past the last observed
  event/censoring time.
  [`ertte_rmst.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
  follows the same convention as
  [`ertte_predict.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_predict.md):
  survival is held flat at its last estimated value beyond that point.
  Because RMST integrates the *entire* curve up to $`\tau`$, this
  flat-tail assumption has a much larger effect on an RMST value than it
  does on a single survival-probability prediction at one time point – a
  long flat tail can inflate the estimated area substantially.
  [`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
  issues a warning whenever any requested `tau` exceeds the observed
  follow-up range, specifically because this matters more here than
  elsewhere in the package.
- For
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md),
  there’s no equivalent hard boundary – the parametric survival curve is
  defined for all $`t > 0`$ – but that comes with its own caveat: far
  beyond the observed data, the estimate is governed entirely by the
  assumed parametric family, with no data left to check whether that
  assumption still holds.

## A worked example

``` r

mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)

nd <- ertte_data[c(1, 50, 100), ]
nd[, c("id", "aucss")]
#>      id    aucss
#> 1     1 1114.089
#> 50   50  148.240
#> 100 100  427.802
```

``` r

ertte_rmst(mod_aft, nd, tau = c(60, 90, 120))
#> # A tibble: 9 × 15
#>      id sex     age weight  dose treatment aucss cmaxss  time event admin_censor
#>   <int> <fct> <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>        <dbl>
#> 1     1 Fema…    27     70   200 Drug      1114.  187.   77.4     1          180
#> 2     1 Fema…    27     70   200 Drug      1114.  187.   77.4     1          180
#> 3     1 Fema…    27     70   200 Drug      1114.  187.   77.4     1          180
#> 4    50 Male     28     58   100 Drug       148.   10.4   0.3     0          180
#> 5    50 Male     28     58   100 Drug       148.   10.4   0.3     0          180
#> 6    50 Male     28     58   100 Drug       148.   10.4   0.3     0          180
#> 7   100 Fema…    26     48   200 Drug       428.   35.1 166.      0          180
#> 8   100 Fema…    26     48   200 Drug       428.   35.1 166.      0          180
#> 9   100 Fema…    26     48   200 Drug       428.   35.1 166.      0          180
#> # ℹ 4 more variables: tau <dbl>, fit_rmst <dbl>, ci_lower <dbl>, ci_upper <dbl>
```

``` r

ertte_rmst(mod_cox, nd, tau = c(60, 90, 120))
#> # A tibble: 9 × 15
#>      id sex     age weight  dose treatment aucss cmaxss  time event admin_censor
#>   <int> <fct> <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>        <dbl>
#> 1     1 Fema…    27     70   200 Drug      1114.  187.   77.4     1          180
#> 2     1 Fema…    27     70   200 Drug      1114.  187.   77.4     1          180
#> 3     1 Fema…    27     70   200 Drug      1114.  187.   77.4     1          180
#> 4    50 Male     28     58   100 Drug       148.   10.4   0.3     0          180
#> 5    50 Male     28     58   100 Drug       148.   10.4   0.3     0          180
#> 6    50 Male     28     58   100 Drug       148.   10.4   0.3     0          180
#> 7   100 Fema…    26     48   200 Drug       428.   35.1 166.      0          180
#> 8   100 Fema…    26     48   200 Drug       428.   35.1 166.      0          180
#> 9   100 Fema…    26     48   200 Drug       428.   35.1 166.      0          180
#> # ℹ 4 more variables: tau <dbl>, fit_rmst <dbl>, ci_lower <dbl>, ci_upper <dbl>
```

The two engines’ point estimates and intervals should generally be in
the same ballpark when both models fit the data reasonably well (as
here, since `ertte_data` was simulated from a Weibull AFT ground truth
that a Cox model should also approximate well) – but they are not
expected to match exactly, and a persistent, large disagreement between
the two is a useful diagnostic that one of the models may not fit the
data well.

``` r

mod_cox_extrap <- ertte_rmst(mod_cox, nd[1, ], tau = 250)
#> Warning: `tau` exceeds the last observed follow-up time (180) for at least one
#> value. RMST beyond that point assumes survival stays flat at its last estimated
#> value, which may understate/overstate the true area.
```

The warning above fires because `250` exceeds the longest observed
follow-up time in `ertte_data`; the returned value still assumes
survival stays flat beyond that point.

## Using RMST in a visual predictive check

`ertte`’s `erplots` integration also supports RMST-based visual
predictive checks (VPCs): passing `tau` through `simulate_args` to
[`erplots::er_vpc_add_simulated()`](https://erplots.djnavarro.net/reference/er_vpc_add_simulated.html)
reduces each simulated replicate to a per-subject restricted survival
time, `min(sim_time, tau)`, then compares its distribution against the
observed data binned by exposure. This reuses the *already-censored*
simulated event times – the same administrative- censoring convention
described above – rather than the model’s raw, uncensored draws, since a
VPC’s purpose is checking whether the model reproduces what you’d
actually observe under the study’s real censoring pattern.

A consequence worth knowing about: a simulated replicate that’s censored
*before* `tau` has a genuinely ambiguous restricted survival time (we
don’t know what would have happened between the censoring time and
`tau`), and is dropped (returned as `NA`) rather than imputed or
reweighted. This is a complete-case convention, applied identically to
the observed and simulated sides of the comparison – and, unlike a naive
complete-case *point estimate*, isn’t obviously the wrong choice here,
since it’s the observed study’s actual censoring pattern driving both
sides equally. A more sophisticated correction
(inverse-probability-of-censoring weighting, or pseudo-observations) was
considered but isn’t offered: erplots’ VPC aggregation is currently an
unweighted mean, so there’s nowhere for an IPCW-style weight to be
applied downstream, and a pseudo-observations approach would need an
expensive leave-one-out recomputation for every simulated replicate. If
this convention drops an uncomfortably large fraction of replicates in a
given analysis, supplying `censor_time` (a per-subject administrative
follow-up horizon) to `simulate_args` changes how much ambiguity arises
in the first place, rather than trying to correct for it after the fact.

## Summary of assumptions and limitations

- **AFT confidence intervals ignore uncertainty in the scale
  parameter**, propagating only the linear predictor’s uncertainty. This
  matches
  [`ertte_predict.ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_predict.md)’s
  existing convention, not a special weakening for RMST specifically.
- **Cox confidence intervals rely on an independent-increments
  approximation** for how coefficient uncertainty accumulates over time,
  which is not an exact result. It has been checked against a bootstrap
  for a small number of covariate profiles and found to track it
  substantially better than the naive alternative, but this is not
  exhaustive validation.
- **Neither engine’s interval is clipped to `[0, tau]`.** Since a
  symmetric Wald interval isn’t automatically bounded the way a
  back-transformed survival-probability interval is,
  `ci_lower`/`ci_upper` can in principle fall outside `[0, tau]` for
  small samples or profiles near the edge of the observed covariate
  range.
  [`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
  does not currently clip this.
- **RMST beyond the observed follow-up range assumes a flat survival
  tail** (Cox engine) or relies entirely on extrapolating the fitted
  parametric family (AFT engine). Treat `tau` values near or beyond the
  longest observed follow-up time with real caution, and note the
  warning
  [`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
  issues for the Cox engine when this happens.
- **RMST depends on the choice of `tau`.** Unlike a hazard ratio, RMST
  is not a single summary of the whole exposure-response relationship –
  it’s specific to the chosen horizon, and different horizons can tell
  different (both valid) stories about the same fitted model.
- **RMST-based VPC simulations drop replicates censored before `tau` as
  ambiguous (complete-case), with no IPCW/pseudo-observation alternative
  currently offered** – see “Using RMST in a visual predictive check”
  above for why, and for the `censor_time` mitigation.

## Further reading

- Royston, P. and Parmar, M.K.B. (2013). Restricted mean survival time:
  an alternative to the hazard ratio for the design and analysis of
  randomized trials with a time-to-event outcome. *BMC Medical Research
  Methodology*, 13, 152.
- [`survival::print.survfit()`](https://rdrr.io/pkg/survival/man/print.survfit.html)’s
  `rmean` argument, and the classic Greenwood-type variance formula it’s
  based on, for the plain (single curve, no covariates) version of this
  idea.
