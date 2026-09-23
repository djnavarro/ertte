# Survival-probability predictions for exposure-response TTE models

Computes fitted survival probabilities `S(t)` and confidence intervals
from a fitted ertte model, for one or more rows of `newdata` at one or
more `time` values.

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

  An ertte model

- ...:

  Passed to methods

- newdata:

  Data frame containing cases to be predicted. Defaults to the data the
  model was fitted to.

- time:

  Numeric vector of times at which to compute survival probabilities

- conf_level:

  Confidence level for the intervals. Defaults to `.95`. Must be a
  single number between 0 and 1 (inclusive); other values error.

## Value

A tibble with one row per combination of `newdata` row and `time`

## Details

`ertte_predict()` is a generic function, with two supplied methods, one
for parametric AFT exposure-response models (i.e., `ertte_aft` classed
objects) and another for Cox proportional hazard exposure-response
models (i.e., objects with class `ertte_coxph`).

For AFT models, it computes the linear predictor and standard error
using [`predict()`](https://rdrr.io/r/stats/predict.html), and then
converts to a survival probability \\S(t)\\:

\$\$S(t) = 1 - F((\log(t) - \mu) / \sigma)\$\$

where \\F(\cdot)\\ denotes the cumulative distribution function for the
AFT base distribution (e.g., Weibull), \\\mu\\ denotes the linear
predictor, and \\\sigma\\ denotes the scale parameter. Confidence
intervals are constructed using Wald intervals on the linear predictor
and back-transformed onto the survival probability scale. Parameter
uncertainty in the scale parameter is not propagated.

For Cox models,
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
is used to compute the survival probability. When `time` exceeds the
last observed follow-up time, the survival function is held constant
(i.e., step-function extrapolation). Confidence intervals are calculated
using the `conf.type = "log"` transform, which specifies Wald intervals
on \\\log(-\log(S))\\, as this is better suited to a probability bounded
in `[0, 1]`. Because of this, it should be noted that the confidence
intervals computed for a Cox model are not directly comparable to those
computed for AFT models.

## Examples

``` r
# predictions for an AFT model
mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
ertte_predict(mod_aft, ertte_data[1:5, ], time = c(30, 60, 90))
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

# predictions for a Cox model
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
