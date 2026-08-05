# Select an AFT distribution by AIC

Fits each candidate AFT distribution to the same formula/data and
selects the best fit by AIC.

## Usage

``` r
ertte_select_distribution(
  formula,
  data,
  candidates = c("exponential", "weibull", "lognormal", "loglogistic")
)
```

## Arguments

- formula:

  Model formula, as for
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)

- data:

  Data set

- candidates:

  Character vector of candidate `dist` values to try. Defaults to all
  four tested/supported distributions.

## Value

A list with two elements: `comparison` (a tibble with one row per
candidate: `dist`, `logLik`, `aic`, `bic`, `converged`, sorted by AIC)
and `model` (the best-fitting `ertte_aft` model, i.e. the one with
lowest AIC).

## Details

Ties (to floating point) are broken by the order `candidates` is given
in, i.e. the first-listed of the tied candidates is returned as `model`.

`candidates` must be a non-empty character vector with no missing
values; each element must also separately name one of the
tested/supported distributions (see
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)'s
`dist` argument). `candidates = character(0)` errors rather than
silently fitting nothing and returning a degenerate
`list(comparison = <0-row tibble>, model = NULL)`.

## Examples

``` r
# (uses ertte_aft() internally for each candidate distribution)
cmp <- ertte_select_distribution(Surv(time, event) ~ aucss, ertte_data)
cmp$comparison
#> # A tibble: 4 × 5
#>   dist        logLik   aic   bic converged
#>   <chr>        <dbl> <dbl> <dbl> <lgl>    
#> 1 weibull     -1207. 2421. 2432. TRUE     
#> 2 loglogistic -1216. 2438. 2449. TRUE     
#> 3 lognormal   -1222. 2451. 2462. TRUE     
#> 4 exponential -1224. 2452. 2460. TRUE     
cmp$model
#> Call:
#> survival::survreg(formula = formula, data = data, dist = dist)
#> 
#> Coefficients:
#>   (Intercept)         aucss 
#>  4.8563375087 -0.0006407913 
#> 
#> Scale= 0.7164032 
#> 
#> Loglik(model)= -1207.3   Loglik(intercept only)= -1263.9
#>  Chisq= 113.25 on 1 degrees of freedom, p= <2e-16 
#> n= 300 
```
