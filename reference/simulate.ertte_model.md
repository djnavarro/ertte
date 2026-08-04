# Simulate from an exposure-response TTE model

[`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
`ertte_model` objects.

## Usage

``` r
# S3 method for class 'ertte_model'
simulate(object, nsim = 100, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  An ertte model object, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)

- nsim:

  Number of simulation replicates

- seed:

  Optional seed. If `NULL` (the default), one is chosen automatically
  and reported via a message (since it determines the actual simulated
  values returned).

- newdata:

  Data frame to simulate from. Defaults to the data the model was fitted
  to. Must contain the original response columns (`time`/`event`, as
  named in the model's `Surv()` call) – see Details.

- ...:

  Unused, present for compatibility with the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) generic

## Value

A tibble with one row per observation per replicate: `dat_id`/`sim_id`,
sampled `coef_*` columns, `sim_time` (the simulated event/censoring
time), and `sim_event` (1 = event, 0 = censored).

## Details

Coefficients are sampled from the asymptotic sampling distribution
implied by `vcov(object)`, and event times are drawn by inverse-CDF
sampling from the fitted AFT distribution at each sampled coefficient
vector (see
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
Details for the log-location-scale representation used). Simulated event
times are capped at each row's *observed* exit time
(`sim_time <- pmin(sim_time_raw, observed_time)`, with `sim_event` set
accordingly) to reproduce the study's observed censoring/follow-up
pattern – a documented simplification, since the true administrative
censoring time for subjects who had an event isn't otherwise available
(see `.ertte_simulate_draws()`).

## Examples

``` r
mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
sim <- simulate(mod, nsim = 20, seed = 1234)
sim
#> # A tibble: 6,000 × 16
#>    dat_id sim_id sim_time sim_event    id sex      age weight  dose treatment
#>     <int>  <int>    <dbl>     <dbl> <int> <fct>  <int>  <dbl> <dbl> <fct>    
#>  1      1      1     77.4         0     1 Female    27     70   200 Drug     
#>  2      2      1     26.8         0     2 Female    27     59   100 Drug     
#>  3      3      1     31.2         1     3 Female    24     65     0 Placebo  
#>  4      4      1     16.8         0     4 Female    29     63     0 Placebo  
#>  5      5      1     17.2         1     5 Male      27     91   200 Drug     
#>  6      6      1     81.3         0     6 Female    18     65     0 Placebo  
#>  7      7      1     46.1         1     7 Male      18     66   200 Drug     
#>  8      8      1     47.9         1     8 Female    20     66   200 Drug     
#>  9      9      1     34.2         1     9 Male      25     62     0 Placebo  
#> 10     10      1    178.          1    10 Male      25     81   100 Drug     
#> # ℹ 5,990 more rows
#> # ℹ 6 more variables: aucss <dbl>, cmaxss <dbl>, time <dbl>, event <dbl>,
#> #   `coef_(Intercept)` <dbl>, coef_aucss <dbl>
```
