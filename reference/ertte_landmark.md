# Landmark event-probability predictions for exposure-response TTE models

Reduces the survival curve for an exposure-response time-to-even model
to a binary landmark event probability at a single fixed time.

## Usage

``` r
ertte_landmark(object, newdata = NULL, landmark_time, conf_level = 0.95, ...)
```

## Arguments

- object:

  An ertte model.

- newdata:

  Data frame containing cases to be predicted. Defaults to the data the
  model was fitted to.

- landmark_time:

  A single positive number.

- conf_level:

  Confidence level for the intervals. Defaults to `.95`.

- ...:

  Passed to
  [`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md).

## Value

A tibble with one row per row of `newdata`. In addition to storing the
original columns from `newdata`, it contains a column storing the
`landmark_time` itself, a `fit_resp` column with the estimated
probability of observing the event at or before the land mark time, and
`ci_lower`, and `ci_upper` columns specifying the confidence interval.

## Details

In some time-to-event analyses it is convenient to reduce a
time-to-event endpoint to a binary landmark response: for every row in
the data, the landmark response describes the modelled probability that
the event occurs at or before a specified `landmark_time` defined by the
analyst. The landmark event probability is simply one minus the survival
probability at that time, and is computed using
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md).
However, unlike the typical usage of
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
in which a vector of times is passed, `landmark_time` must be a single
fixed value – a landmark is by definition evaluated at one time.

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
