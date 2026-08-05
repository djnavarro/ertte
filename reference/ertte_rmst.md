# Restricted mean survival time predictions for exposure-response TTE models

Restricted mean survival time predictions for exposure-response TTE
models

## Usage

``` r
ertte_rmst(object, newdata = NULL, tau, conf_level = 0.95, ...)

# S3 method for class 'ertte_aft'
ertte_rmst(object, newdata = NULL, tau, conf_level = 0.95, ...)

# S3 method for class 'ertte_coxph'
ertte_rmst(object, newdata = NULL, tau, conf_level = 0.95, ...)
```

## Arguments

- object:

  An ertte model, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md).

- newdata:

  Data frame containing cases to be predicted. Defaults to the data the
  model was fitted to.

- tau:

  Numeric vector of restriction horizons at which to compute
  `RMST(tau) = integral of S(t) from 0 to tau`.

- conf_level:

  Confidence level for the intervals.

- ...:

  Passed to methods.

## Value

A tibble with one row per combination of `newdata` row and `tau`, plus
`fit_rmst`, `ci_lower`, and `ci_upper`.

## Details

`ertte_rmst()` reduces a TTE endpoint to a scalar exposure-response
value – restricted mean survival time, the area under the survival curve
up to a fixed horizon `tau` – the other scalar reduction the package's
design issue mentions alongside landmark-binary (see
[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)).
Unlike
[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md),
`ertte_rmst()` is a generic (not a thin wrapper around
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)):
computing an area under the curve genuinely needs engine-specific logic
(the whole curve, not a single time point), so there are per-engine
methods – see `ertte_rmst.ertte_aft()` and `ertte_rmst.ertte_coxph()`.

Confidence intervals are symmetric Wald intervals on the RMST scale
(`fit_rmst +/- z * se_rmst`) for both engines. This is not automatically
bounded to `[0, tau]` the way
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)'s
survival-probability intervals are bounded to `[0, 1]` by construction
(their back-transform through a CDF keeps them there) – an unclipped
Wald interval on RMST can, in principle, dip below 0 or exceed `tau` for
small samples or near-boundary cases. `conf_level` must be a single
number between 0 and 1 (inclusive); other values error rather than
silently producing a reversed or `NaN` interval.

A zero-row `newdata` returns a zero-row tibble with the expected columns
for both engines, rather than erroring (see issue \#10 for why this
needed an explicit guard on the `ertte_coxph` side).

The `ertte_aft` method computes `fit_rmst` by numerically integrating
the closed-form survival function \`S(t\|x) = 1 - F((log(t)

- mu) /
  scale)`from 0 to`tau`via`stats::integrate()`, where `mu`(and its standard error) comes from`predict(object,
  newdata, type = "linear", se.fit =
  TRUE)`, matching [ertte_predict.ertte_aft()]. The standard error is an analytic delta method that differentiates under the integral sign: `d/dmu
  RMST(tau\|x) = integral of dbase(z) / scale from 0 to
  tau`, where `dbase`is the base distribution's density (see`.ertte_dist_info()`) -- propagating only `Var(mu)`, not `Var(scale)`, the same simplification `ertte_predict.ertte_aft()\`
  already makes for its own confidence intervals.

The `ertte_coxph` method delegates to
`survival::survfit(object, newdata, conf.int = conf_level, se.fit = TRUE)`,
the same call
[`ertte_predict.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
uses. Because the fitted baseline hazard (and therefore every
covariate-adjusted survival curve) is a right-continuous step function,
`fit_rmst` is an *exact* finite sum of rectangle areas between
consecutive jump times up to `tau` – not a numerical-quadrature
approximation. `se_rmst` comes from a delta method that reuses the same
rectangle/tail-weighted construction as `survival:::survmean()`'s
classic Greenwood-based RMST variance, but with the variance-increment
term replaced by increments of the profile-specific `std.err(t)^2`
returned by `survfit()` – `survmean()`'s own Greenwood term is based on
population-level risk sets shared across every covariate profile, which
understates uncertainty for a profile far from the mean covariate values
(see `.ertte_rmst_pfun_delta()`'s source comments for the derivation and
the bootstrap cross-check that motivated this).

If any value of `tau` exceeds the last observed follow-up time across
the whole fitted cohort, `ertte_rmst()` warns: RMST integrates the
*entire* curve up to `tau`, so silently assuming survival stays flat
beyond the observed range (the same extrapolation convention
[`ertte_predict.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
uses for a single time point) has a larger effect on an area than on a
point-in-time prediction.

`conf_level = 0`/`1`, documented (see `.ertte_check_conf_level()`) as
legitimate degenerate endpoints, are supported here directly, since the
delta-method interval is built from `z_scale`
([`qnorm()`](https://rdrr.io/r/stats/Normal.html)-derived, `0` or `Inf`
at these boundaries) rather than `survfit()`'s own `conf.int` machinery
– the latter is only used to request `$surv`/ `$std.err`, which don't
depend on the requested confidence level (see issue \#11); `survfit()`
is always called with a fixed, valid placeholder value internally.

A zero-row `newdata` returns a zero-row tibble with the expected columns
rather than erroring:
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
itself rejects an entirely-missing `newdata` with a cryptic "all rows of
newdata have missing values" error (see issue \#10).

## Examples

``` r
mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
ertte_rmst(mod, ertte_data[1:5, ], tau = c(60, 90))
#> # A tibble: 10 × 15
#>       id sex      age weight  dose treatment aucss cmaxss  time event
#>    <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>
#>  1     1 Female    27     70   200 Drug      1114.  187.   77.4     1
#>  2     1 Female    27     70   200 Drug      1114.  187.   77.4     1
#>  3     2 Female    27     59   100 Drug       561.   49.1  26.8     0
#>  4     2 Female    27     59   100 Drug       561.   49.1  26.8     0
#>  5     3 Female    24     65     0 Placebo      0     0   180       0
#>  6     3 Female    24     65     0 Placebo      0     0   180       0
#>  7     4 Female    29     63     0 Placebo      0     0    16.8     0
#>  8     4 Female    29     63     0 Placebo      0     0    16.8     0
#>  9     5 Male      27     91   200 Drug      1416.  143.   33.9     1
#> 10     5 Male      27     91   200 Drug      1416.  143.   33.9     1
#> # ℹ 5 more variables: admin_censor <dbl>, tau <dbl>, fit_rmst <dbl>,
#> #   ci_lower <dbl>, ci_upper <dbl>

mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
ertte_rmst(mod_cox, ertte_data[1:5, ], tau = c(60, 90))
#> # A tibble: 10 × 15
#>       id sex      age weight  dose treatment aucss cmaxss  time event
#>    <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>
#>  1     1 Female    27     70   200 Drug      1114.  187.   77.4     1
#>  2     1 Female    27     70   200 Drug      1114.  187.   77.4     1
#>  3     2 Female    27     59   100 Drug       561.   49.1  26.8     0
#>  4     2 Female    27     59   100 Drug       561.   49.1  26.8     0
#>  5     3 Female    24     65     0 Placebo      0     0   180       0
#>  6     3 Female    24     65     0 Placebo      0     0   180       0
#>  7     4 Female    29     63     0 Placebo      0     0    16.8     0
#>  8     4 Female    29     63     0 Placebo      0     0    16.8     0
#>  9     5 Male      27     91   200 Drug      1416.  143.   33.9     1
#> 10     5 Male      27     91   200 Drug      1416.  143.   33.9     1
#> # ℹ 5 more variables: admin_censor <dbl>, tau <dbl>, fit_rmst <dbl>,
#> #   ci_lower <dbl>, ci_upper <dbl>
```
