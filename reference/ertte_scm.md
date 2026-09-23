# Stepwise covariate modelling for exposure-response TTE models

Iteratively adds (`ertte_scm_forward()`) or removes
(`ertte_scm_backward()`) covariate terms from an ertte model one at a
time, selecting each step's term by a likelihood-ratio p-value, AIC, or
BIC, and logging the search in a history object retrievable via
`ertte_scm_history()`.

## Usage

``` r
ertte_scm_forward(
  mod,
  candidates,
  threshold = 0.01,
  criterion = "p-value",
  seed = NULL
)

ertte_scm_backward(
  mod,
  candidates,
  threshold = 0.001,
  criterion = "p-value",
  seed = NULL
)

ertte_scm_history(mod)
```

## Arguments

- mod:

  An ertte model object

- candidates:

  Character vector with list of candidate terms

- threshold:

  Threshold to test against. Used only when `criterion = "p-value"` (the
  default); ignored otherwise. Defaults to `0.01` for
  `ertte_scm_forward()` and `0.001` for `ertte_scm_backward()`.

- criterion:

  Model selection criterion. One of `"p-value"` (default), `"aic"`, or
  `"bic"`.

- seed:

  Optional seed to control the order candidate terms are tested within a
  step. If `NULL` (the default), one is chosen automatically (see
  "Candidate test order and `seed`" below for when this actually
  matters).

## Value

For `ertte_scm_forward()` and `ertte_scm_backward()`, the updated ertte
model is returned, with the SCM history log updated internally. For
`ertte_scm_history()`, a data frame is returned containing the SCM
history log.

## Details

Three model selection criteria are available via the `criterion`
argument:

- `"p-value"` (default): Models are compared with a likelihood-ratio
  Chi-squared test. A term is added if its likelihood-ratio p-value
  falls below `threshold` (forward) or removed if its p-value exceeds
  `threshold` (backward). When multiple candidates satisfy the threshold
  within a step, the one with the most extreme p-value is chosen.

- `"aic"`: A term is added (forward) or removed (backward) if doing so
  strictly decreases AIC relative to the current model. When multiple
  candidates improve AIC, the one yielding the lowest AIC is chosen.

- `"bic"`: Same as `"aic"`, but using BIC as the criterion.

When `criterion` is `"aic"` or `"bic"`, the `threshold` argument has no
effect and is ignored, and `term_p_value` is left `NA` in the history
for every candidate tested that step (the likelihood-ratio test isn't
computed, since it plays no role in selection). The `model_aic` and
`model_bic` columns are always recorded regardless of which criterion
drove selection, and the history's `criterion` column records which one
was used for each forward/backward step.

## Handling problem candidates

A candidate that can't be fit for a given step (e.g. it's aliased with a
term already in the model, references a variable missing from the data,
or fails to converge) is skipped for that step with a warning, rather
than aborting the whole search or being silently selected.

## Candidate test order and `seed`

`seed` exists as a safety measure against run-to-run variation in the
order candidate terms are tested within a step
([`sample()`](https://rdrr.io/r/base/sample.html), shuffled before
testing one at a time). Model fitting itself
([`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html))
is deterministic given a starting formula, so `seed` only matters in the
(essentially measure-zero) case of an exact p-value tie between
competing candidates within a step – see the companion `erglm` package's
equivalent documentation for the full rationale, which applies unchanged
here.

## Examples

``` r
# Forward addition by p-value selection
mod0 <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"))
ertte_scm_history(mod1)
#> # A tibble: 4 × 12
#>   iteration attempt step       criterion action term_tested model_tested        
#>       <int>   <int> <chr>      <chr>     <chr>  <chr>       <chr>               
#> 1         0       0 base model NA        NA     NA          Surv(time, event) ~…
#> 2         1       1 forward    p-value   add    ~sex        Surv(time, event) ~…
#> 3         1       2 forward    p-value   add    ~dose       Surv(time, event) ~…
#> 4         2       3 forward    p-value   add    ~dose       Surv(time, event) ~…
#> # ℹ 5 more variables: model_converged <lgl>, term_p_value <dbl>,
#> #   model_aic <dbl>, model_bic <dbl>, model_updated <int>

# Backward elimination by p-value selection
mod2 <- ertte_aft(Surv(time, event) ~ aucss + sex + dose, ertte_data)
mod3 <- ertte_scm_backward(mod2, candidates = c("sex", "dose"))
ertte_scm_history(mod3)
#> # A tibble: 4 × 12
#>   iteration attempt step       criterion action term_tested model_tested        
#>       <int>   <int> <chr>      <chr>     <chr>  <chr>       <chr>               
#> 1         0       0 base model NA        NA     NA          Surv(time, event) ~…
#> 2         1       1 backward   p-value   remove ~dose       Surv(time, event) ~…
#> 3         1       2 backward   p-value   remove ~sex        Surv(time, event) ~…
#> 4         2       3 backward   p-value   remove ~sex        Surv(time, event) ~…
#> # ℹ 5 more variables: model_converged <lgl>, term_p_value <dbl>,
#> #   model_aic <dbl>, model_bic <dbl>, model_updated <int>

# AIC-based forward addition/backward elimination instead of p-value
mod4 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"), criterion = "aic")
mod5 <- ertte_scm_backward(mod4, candidates = c("sex", "dose"), criterion = "bic")
ertte_scm_history(mod5)
#> # A tibble: 5 × 12
#>   iteration attempt step       criterion action term_tested model_tested        
#>       <int>   <int> <chr>      <chr>     <chr>  <chr>       <chr>               
#> 1         0       0 base model NA        NA     NA          Surv(time, event) ~…
#> 2         1       1 forward    aic       add    ~dose       Surv(time, event) ~…
#> 3         1       2 forward    aic       add    ~sex        Surv(time, event) ~…
#> 4         2       3 forward    aic       add    ~dose       Surv(time, event) ~…
#> 5         3       4 backward   bic       remove ~sex        Surv(time, event) ~…
#> # ℹ 5 more variables: model_converged <lgl>, term_p_value <dbl>,
#> #   model_aic <dbl>, model_bic <dbl>, model_updated <int>
```
