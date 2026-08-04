# Simulate from an exposure-response TTE model

[`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
`ertte_model` objects. Works for both `ertte_aft` and `ertte_coxph` fits
via a single shared method – there's no separate
`simulate.ertte_coxph()` – since the engine-specific simulation
mechanics are dispatched internally by `.ertte_simulate_draws()` (see
Details).

## Usage

``` r
# S3 method for class 'ertte_model'
simulate(
  object,
  nsim = 100,
  seed = NULL,
  newdata = NULL,
  censor_time = NULL,
  ...
)
```

## Arguments

- object:

  An ertte model object, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)

- nsim:

  Number of simulation replicates

- seed:

  Optional seed. If `NULL` (the default), one is chosen automatically
  and reported via a message (since it determines the actual simulated
  values returned).

- newdata:

  Data frame to simulate from. Defaults to the data the model was fitted
  to. Must contain the original response columns (`time`/`event`, as
  named in the model's
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) call) – see
  Details.

- censor_time:

  Optional administrative/maximum-follow-up time(s) to cap simulated
  event times at, applied uniformly to every row regardless of whether
  that row observed an event. Either `NULL` (the default – see Details
  for the fallback behaviour), a single number (recycled across all rows
  of `newdata`), or a numeric vector of length `nrow(newdata)` giving
  each row's own administrative follow-up time.

- ...:

  Unused, present for compatibility with the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) generic

## Value

A tibble with one row per observation per replicate: `dat_id`/`sim_id`,
sampled `coef_*` columns, `sim_time` (the simulated event/censoring
time), and `sim_event` (1 = event, 0 = censored).

## Details

Coefficients are sampled from the asymptotic sampling distribution
implied by `vcov(object)`. Event times are then drawn by inverse-CDF
sampling, via the internal `.ertte_simulate_draws()` S3 generic, whose
per-engine methods differ in exactly how: for `ertte_aft` fits, directly
from the fitted log-location-scale AFT distribution (see
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
Details); for `ertte_coxph` fits, by inverting the fitted baseline
cumulative hazard
([`survival::basehaz()`](https://rdrr.io/pkg/survival/man/basehaz.html),
held fixed regardless of the sampled coefficient draw – the same
simplification
[`ertte_fun.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_fun.md)
makes for a user-supplied `param`).

The resulting *raw* (uncensored) simulated event time is then censored
by `censor_time` if supplied
(`sim_time <- pmin(sim_time_raw, censor_time)`, with `sim_event` set
accordingly) – this is the accurate case, whenever a genuine per-row (or
study-wide constant) administrative follow-up time is known, since it
caps every row (event or censored) against its true censoring horizon.

Absent a supplied `censor_time` (the default, `NULL`), row's own
observed `event` status determines the fallback: rows that were
*censored* in `newdata` have their observed exit time used as the cap
(`sim_time <- pmin(sim_time_raw, observed_time)`), since that observed
exit time genuinely is when censoring happened – an exact match, not an
approximation. Rows that had an observed *event*, however, are left
**uncensored** in the simulation: their observed exit time is when the
event actually happened, not their administrative censoring horizon
(which was necessarily later, and typically isn't recorded once an event
has occurred) – capping simulated draws there would leak the observed
event day into the simulation and bias a simulated-vs-observed
comparison (e.g. a visual predictive check) toward looking more similar
than the fitted model actually implies. This remains an approximation
for event rows (no censoring is applied at all, absent better
information), but avoids that specific bias – see
`.ertte_apply_admin_censoring()`.

## Examples

``` r
mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
sim <- simulate(mod, nsim = 20, seed = 1234)
sim
#> # A tibble: 6,000 × 17
#>    dat_id sim_id sim_time sim_event    id sex      age weight  dose treatment
#>     <int>  <int>    <dbl>     <dbl> <int> <fct>  <int>  <dbl> <dbl> <fct>    
#>  1      1      1    122.          1     1 Female    27     70   200 Drug     
#>  2      2      1     26.8         0     2 Female    27     59   100 Drug     
#>  3      3      1     31.2         1     3 Female    24     65     0 Placebo  
#>  4      4      1     16.8         0     4 Female    29     63     0 Placebo  
#>  5      5      1     17.2         1     5 Male      27     91   200 Drug     
#>  6      6      1    216.          1     6 Female    18     65     0 Placebo  
#>  7      7      1     46.1         1     7 Male      18     66   200 Drug     
#>  8      8      1     47.9         1     8 Female    20     66   200 Drug     
#>  9      9      1     34.2         1     9 Male      25     62     0 Placebo  
#> 10     10      1    178.          1    10 Male      25     81   100 Drug     
#> # ℹ 5,990 more rows
#> # ℹ 7 more variables: aucss <dbl>, cmaxss <dbl>, time <dbl>, event <dbl>,
#> #   admin_censor <dbl>, `coef_(Intercept)` <dbl>, coef_aucss <dbl>

mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
sim_cox <- simulate(mod_cox, nsim = 20, seed = 1234)
sim_cox
#> # A tibble: 6,000 × 16
#>    dat_id sim_id sim_time sim_event    id sex      age weight  dose treatment
#>     <int>  <int>    <dbl>     <dbl> <int> <fct>  <int>  <dbl> <dbl> <fct>    
#>  1      1      1     46.4         1     1 Female    27     70   200 Drug     
#>  2      2      1     26.8         0     2 Female    27     59   100 Drug     
#>  3      3      1    135.          1     3 Female    24     65     0 Placebo  
#>  4      4      1     16.8         0     4 Female    29     63     0 Placebo  
#>  5      5      1     56           1     5 Male      27     91   200 Drug     
#>  6      6      1     88.4         1     6 Female    18     65     0 Placebo  
#>  7      7      1     43.6         1     7 Male      18     66   200 Drug     
#>  8      8      1     74.1         1     8 Female    20     66   200 Drug     
#>  9      9      1    161.          1     9 Male      25     62     0 Placebo  
#> 10     10      1     43.2         1    10 Male      25     81   100 Drug     
#> # ℹ 5,990 more rows
#> # ℹ 6 more variables: aucss <dbl>, cmaxss <dbl>, time <dbl>, event <dbl>,
#> #   admin_censor <dbl>, coef_aucss <dbl>

# a genuine per-row administrative censoring time -- ertte_data's
# `admin_censor` column records the fixed 180-day study cutoff used
# to generate it, known regardless of whether a subject had an event
sim_admin <- simulate(mod, nsim = 20, seed = 1234, censor_time = ertte_data$admin_censor)
sim_admin
#> # A tibble: 6,000 × 17
#>    dat_id sim_id sim_time sim_event    id sex      age weight  dose treatment
#>     <int>  <int>    <dbl>     <dbl> <int> <fct>  <int>  <dbl> <dbl> <fct>    
#>  1      1      1    122.          1     1 Female    27     70   200 Drug     
#>  2      2      1     61.9         1     2 Female    27     59   100 Drug     
#>  3      3      1     31.2         1     3 Female    24     65     0 Placebo  
#>  4      4      1    100.          1     4 Female    29     63     0 Placebo  
#>  5      5      1     17.2         1     5 Male      27     91   200 Drug     
#>  6      6      1    180           0     6 Female    18     65     0 Placebo  
#>  7      7      1     46.1         1     7 Male      18     66   200 Drug     
#>  8      8      1     47.9         1     8 Female    20     66   200 Drug     
#>  9      9      1     34.2         1     9 Male      25     62     0 Placebo  
#> 10     10      1    178.          1    10 Male      25     81   100 Drug     
#> # ℹ 5,990 more rows
#> # ℹ 7 more variables: aucss <dbl>, cmaxss <dbl>, time <dbl>, event <dbl>,
#> #   admin_censor <dbl>, `coef_(Intercept)` <dbl>, coef_aucss <dbl>
```
