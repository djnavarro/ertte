# Survival-probability predictions for exposure-response TTE models

Survival-probability predictions for exposure-response TTE models

## Usage

``` r
ertte_predict(object, newdata = NULL, time, conf_level = 0.95)
```

## Arguments

- object:

  An ertte model, as returned by
  [`ertte_model()`](https://ertte.djnavarro.net/reference/ertte_model.md)

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

Computes the linear predictor (and its standard error) via
`predict(object, newdata, type = "linear", se.fit = TRUE)`, then
converts to a survival probability
`S(t) = 1 - F((log(t) - mu) / scale)`, where `F` is the base
distribution's CDF implied by `object`'s `dist` (see
[`ertte_model()`](https://ertte.djnavarro.net/reference/ertte_model.md)
Details). Confidence intervals are Wald intervals on `mu` (a
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) z-score times the
standard error), back-transformed the same way – parameter uncertainty
in `scale` is not propagated, matching the level of approximation used
throughout this package (e.g. `erglm_predict()`'s equivalent in the
companion `erglm` package). `conf_level` must be a single number between
0 and 1 (inclusive); other values error rather than silently producing a
reversed or `NaN` interval.

## Examples

``` r
mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
ertte_predict(mod, ertte_data[1:5, ], time = c(30, 60, 90))
#> # A tibble: 15 × 13
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
#> # ℹ 3 more variables: fit_survival <dbl>, ci_lower <dbl>, ci_upper <dbl>
```
