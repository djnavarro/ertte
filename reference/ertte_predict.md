# Survival-probability predictions for exposure-response TTE models

Survival-probability predictions for exposure-response TTE models

## Usage

``` r
ertte_predict(object, ...)

# S3 method for class 'ertte_aft'
ertte_predict(object, newdata = NULL, time, conf_level = 0.95, ...)

# S3 method for class 'ertte_coxph'
ertte_predict(object, newdata = NULL, time, conf_level = 0.95, ...)
```

## Arguments

- object:

  An ertte model, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  (not yet implemented)

- ...:

  Passed to methods

- newdata:

  Data frame containing cases to be predicted. Defaults to the data the
  model was fitted to.

- time:

  Numeric vector of times at which to compute survival probabilities

- conf_level:

  Confidence level for the intervals

## Value

A tibble with one row per combination of `newdata` row and `time`

## Details

`ertte_predict()` is a generic, with methods for each supported engine –
see `ertte_predict.ertte_aft()`.

The `ertte_aft` method computes the linear predictor (and its standard
error) via `predict(object, newdata, type = "linear", se.fit = TRUE)`,
then converts to a survival probability
`S(t) = 1 - F((log(t) - mu) / scale)`, where `F` is the base
distribution's CDF implied by `object`'s `dist` (see
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
Details). Confidence intervals are Wald intervals on `mu` (a
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) z-score times the
standard error), back-transformed the same way – parameter uncertainty
in `scale` is not propagated, matching the level of approximation used
throughout this package (e.g. `erglm_predict()`'s equivalent in the
companion `erglm` package). `conf_level` must be a single number between
0 and 1 (inclusive); other values error rather than silently producing a
reversed or `NaN` interval.

The `ertte_coxph` method delegates to
`survival::survfit(object, newdata, conf.int = conf_level)`, which
computes a per-row survival curve
`S(t | x) = S0(t)^exp(lp(x) - lp(xbar))` from the fitted baseline hazard
(Breslow or Efron, matching `object$method`) and the linear predictor,
then evaluates it at `time` via `summary(..., extend = TRUE)` –
`extend = TRUE` allows `time` to exceed the last observed follow-up
time, holding survival constant beyond it (the usual step-function
extrapolation) rather than erroring. Confidence intervals come from
`survfit()`'s own `conf.type = "log"` transform (Wald on
`log(-log(S))`), which is better suited to a probability bounded in
`[0, 1]` than the plain Wald interval `ertte_predict.ertte_aft()` uses
on the linear predictor – the two methods' intervals are not directly
comparable as a result, which is expected given the different model
structures.

## Examples

``` r
mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
ertte_predict(mod, ertte_data[1:5, ], time = c(30, 60, 90))
#> # A tibble: 15 × 14
#>       id sex      age weight  dose treatment aucss cmaxss  time event
#>    <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>
#>  1     1 Female    27     70   200 Drug      1114.  187.     30     1
#>  2     1 Female    27     70   200 Drug      1114.  187.     60     1
#>  3     1 Female    27     70   200 Drug      1114.  187.     90     1
#>  4     2 Female    27     59   100 Drug       561.   49.1    30     0
#>  5     2 Female    27     59   100 Drug       561.   49.1    60     0
#>  6     2 Female    27     59   100 Drug       561.   49.1    90     0
#>  7     3 Female    24     65     0 Placebo      0     0      30     0
#>  8     3 Female    24     65     0 Placebo      0     0      60     0
#>  9     3 Female    24     65     0 Placebo      0     0      90     0
#> 10     4 Female    29     63     0 Placebo      0     0      30     0
#> 11     4 Female    29     63     0 Placebo      0     0      60     0
#> 12     4 Female    29     63     0 Placebo      0     0      90     0
#> 13     5 Male      27     91   200 Drug      1416.  143.     30     1
#> 14     5 Male      27     91   200 Drug      1416.  143.     60     1
#> 15     5 Male      27     91   200 Drug      1416.  143.     90     1
#> # ℹ 4 more variables: admin_censor <dbl>, fit_survival <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>

mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
ertte_predict(mod_cox, ertte_data[1:5, ], time = c(30, 60, 90))
#> # A tibble: 15 × 14
#>       id sex      age weight  dose treatment aucss cmaxss  time event
#>    <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>
#>  1     1 Female    27     70   200 Drug      1114.  187.     30     1
#>  2     1 Female    27     70   200 Drug      1114.  187.     60     1
#>  3     1 Female    27     70   200 Drug      1114.  187.     90     1
#>  4     2 Female    27     59   100 Drug       561.   49.1    30     0
#>  5     2 Female    27     59   100 Drug       561.   49.1    60     0
#>  6     2 Female    27     59   100 Drug       561.   49.1    90     0
#>  7     3 Female    24     65     0 Placebo      0     0      30     0
#>  8     3 Female    24     65     0 Placebo      0     0      60     0
#>  9     3 Female    24     65     0 Placebo      0     0      90     0
#> 10     4 Female    29     63     0 Placebo      0     0      30     0
#> 11     4 Female    29     63     0 Placebo      0     0      60     0
#> 12     4 Female    29     63     0 Placebo      0     0      90     0
#> 13     5 Male      27     91   200 Drug      1416.  143.     30     1
#> 14     5 Male      27     91   200 Drug      1416.  143.     60     1
#> 15     5 Male      27     91   200 Drug      1416.  143.     90     1
#> # ℹ 4 more variables: admin_censor <dbl>, fit_survival <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>
```
