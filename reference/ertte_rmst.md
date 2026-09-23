# Restricted mean survival time predictions for exposure-response TTE models

Computes restricted mean survival time (RMST), the area under a fitted
ertte model's survival curve up to one or more fixed horizons `tau`,
with confidence intervals.

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

  Numeric vector of restriction horizons at which to compute RSMT.

- conf_level:

  Confidence level for the intervals. Defaults to `.95`.

- ...:

  Passed to methods.

## Value

A tibble with one row per combination of `newdata` row and `tau`, plus
`fit_rmst`, `ci_lower`, and `ci_upper`.

## Details

`ertte_rmst()` is used to reduce a time-to-event endpoint to a scalar
exposure-response value, in this case the restricted mean survival time
(RMST). The RMST for some horizon time `tau` is defined as the area
under the survival function calculated up to the horizon:

\$\$\mbox{RMST}(\tau) = \int_0^\tau S(t) \\ dt\$\$

Because this calculation is performed differently for parametric AFT
models and the Cox semiparametric proportional hazard model,
`ertte_rmst()` is a generic function with two defined methods, one for
each kind of model supported by the ertte package.

For both model types, confidence intervals are symmetric Wald intervals
on the RMST scale. This is not automatically bounded to `[0, tau]` the
way that the survival probability intervals for
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
are bounded to `[0, 1]`. An unclipped Wald interval on RMST can, in
principle, dip below 0 or exceed `tau` for small samples or
near-boundary cases.

If any value of `tau` exceeds the last observed follow-up time across
the whole fitted cohort, `ertte_rmst()` produces a warning that informs
the user that the assumption that survival stays flat beyond the
observed range may be unreliable as this assumption has a larger effect
on an area than on a point-in-time prediction.

## AFT method

The `ertte_aft` method computes `fit_rmst` by numerically integrating
the closed-form survival function via
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html). The
integration is performed on the log-time scale rather than directly
integrating over time, as adaptive quadrature can fail silently on the
raw time scale. Even with this, however, the integration is not always
reliable for extreme values of `tau`, and so `ertte_rmst()` produces
warnings when `tau` is especially large.

## Cox PH method

The `ertte_coxph` method delegates the work to
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
Because the fitted baseline hazard (and therefore every
covariate-adjusted survival curve) is a right-continuous step function,
`fit_rmst` is an *exact* finite sum of rectangle areas between
consecutive jump times up to `tau`: it is not a numerical-quadrature
approximation. The standard error column `se_rmst` comes from a delta
method that reuses the same rectangle/tail-weighted construction as the
Greenwood-based RMST variance method in `survival:::survmean()`, but
with the variance-increment term replaced by increments of the
profile-specific `std.err(t)^2` returned by `survfit()`. In contrast,
`survmean()`'s own Greenwood term is based on population-level risk sets
shared across every covariate profile, which understates uncertainty for
a profile far from the mean covariate values.

## Examples

``` r
# RMST for a parametric AFT model
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

# RMST for a semiparametric Cox PH model
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
