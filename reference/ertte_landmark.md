# Landmark event-probability predictions for exposure-response TTE models

Landmark event-probability predictions for exposure-response TTE models

## Usage

``` r
ertte_landmark(object, newdata = NULL, landmark_time, conf_level = 0.95, ...)
```

## Arguments

- object:

  An ertte model, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md).

- newdata:

  Data frame containing cases to be predicted. Defaults to the data the
  model was fitted to.

- landmark_time:

  A single, strictly positive number: the fixed time `t*` at which to
  compute `P(event by t*)`.

- conf_level:

  Confidence level for the intervals.

- ...:

  Passed to
  [`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md).

## Value

A tibble with one row per row of `newdata`, plus `landmark_time`,
`fit_resp` (the estimated `P(event by t*)`), `ci_lower`, and `ci_upper`.

## Details

Reduces a TTE endpoint to a binary landmark response – "did the event
happen by a fixed time t\*" – turning it into an ordinary scalar
exposure-response value that erplots' existing `er_plot()`/`er_vpc()`
grammars can visualise with no new plotting code (see the package's
design issue, Workstream B1:
<https://github.com/djnavarro/ertte/issues/1>). \`P(event by t\*) = 1

- S(t\*)`, computed by calling [ertte_predict()] at `time =
  landmark_time`and transforming its survival-probability output. Since that's a decreasing monotonic transform, the confidence interval bounds swap (the upper bound on survival becomes the lower bound on event probability, and vice versa) but need no recomputation of their own: whatever validity`ertte_predict()`'s interval has for a given engine -- a Wald interval on the AFT method's linear predictor, or `survival::survfit()`'s own `conf.type
  = "log"\` interval for the Cox PH method – carries through unchanged.

`ertte_landmark()` is a single function, not a generic – unlike
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)/[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md),
it needs no engine-specific logic of its own: it delegates entirely to
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md),
which already dispatches on the `ertte_aft`/`ertte_coxph` subclass. This
also means all of
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)'s
existing edge-case handling (e.g. the all-censored-Cox guard,
single-stratum `NA` propagation) is inherited unchanged.

Unlike
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)'s
`time` argument (a vector, evaluated at potentially many times per row),
`landmark_time` must be a single fixed value – a landmark is by
definition evaluated at one time.

Restricted mean survival time (RMST), the other scalar E-R reduction the
design issue mentions, is not implemented here (see AGENTS.md's "Planned
work").

## Examples

``` r
mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
ertte_landmark(mod, ertte_data[1:5, ], landmark_time = 180)
#> # A tibble: 5 × 14
#>      id sex      age weight  dose treatment aucss cmaxss event admin_censor
#>   <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl>        <dbl>
#> 1     1 Female    27     70   200 Drug      1114.  187.      1          180
#> 2     2 Female    27     59   100 Drug       561.   49.1     0          180
#> 3     3 Female    24     65     0 Placebo      0     0       0          180
#> 4     4 Female    29     63     0 Placebo      0     0       0          180
#> 5     5 Male      27     91   200 Drug      1416.  143.      1          180
#> # ℹ 4 more variables: landmark_time <dbl>, fit_resp <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>

mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
ertte_landmark(mod_cox, ertte_data[1:5, ], landmark_time = 180)
#> # A tibble: 5 × 14
#>      id sex      age weight  dose treatment aucss cmaxss event admin_censor
#>   <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl>        <dbl>
#> 1     1 Female    27     70   200 Drug      1114.  187.      1          180
#> 2     2 Female    27     59   100 Drug       561.   49.1     0          180
#> 3     3 Female    24     65     0 Placebo      0     0       0          180
#> 4     4 Female    29     63     0 Placebo      0     0       0          180
#> 5     5 Male      27     91   200 Drug      1416.  143.      1          180
#> # ℹ 4 more variables: landmark_time <dbl>, fit_resp <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>
```
