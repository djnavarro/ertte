# Sample simulated data for exposure-response time-to-event models

Sample simulated data for exposure-response time-to-event models

## Usage

``` r
ertte_data
```

## Format

A data frame with columns:

- id:

  Identifier

- sex:

  Sex

- age:

  Age

- weight:

  Weight

- dose:

  Nominal dose, units not specified

- treatment:

  Treatment

- aucss:

  AUCss

- cmaxss:

  Cmax,ss

- time:

  Time to event or right-censoring, in days

- event:

  Event indicator (1 = event observed, 0 = right-censored)

- admin_censor:

  Administrative censoring cutoff (days) – fixed at 180 for every
  subject, regardless of whether that subject had an event. Demonstrates
  the `censor_time` argument to
  [`simulate.ertte_model()`](https://ertte.djnavarro.net/reference/simulate.ertte_model.md),
  since (unlike `time`) it's knowable for event rows too.

## Details

This simulated dataset is entirely synthetic. `time`/`event` were
generated from a Weibull AFT model with an exposure (`aucss`) and sex
effect, administrative censoring at 180 days (retained as
`admin_censor`), and a little independent random (dropout) censoring
(not retained, since – like most real TTE data – a subject's dropout
time isn't observable once an event has occurred). You can find the data
generating code in the package source code.

## Examples

``` r
ertte_data
#> # A tibble: 300 × 11
#>       id sex      age weight  dose treatment aucss cmaxss  time event
#>    <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>
#>  1     1 Female    27     70   200 Drug      1114.  187.   77.4     1
#>  2     2 Female    27     59   100 Drug       561.   49.1  26.8     0
#>  3     3 Female    24     65     0 Placebo      0     0   180       0
#>  4     4 Female    29     63     0 Placebo      0     0    16.8     0
#>  5     5 Male      27     91   200 Drug      1416.  143.   33.9     1
#>  6     6 Female    18     65     0 Placebo      0     0    81.3     1
#>  7     7 Male      18     66   200 Drug       748.   42.1  55.5     1
#>  8     8 Female    20     66   200 Drug       344.   35.6 123.      1
#>  9     9 Male      25     62     0 Placebo      0     0   165.      1
#> 10    10 Male      25     81   100 Drug       313.   13.7 180       0
#> # ℹ 290 more rows
#> # ℹ 1 more variable: admin_censor <dbl>
```
