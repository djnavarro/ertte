# Stepwise covariate modelling for exposure-response TTE models

Stepwise covariate modelling for exposure-response TTE models

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
  default); ignored otherwise.

- criterion:

  Model selection criterion. One of `"p-value"` (default), `"aic"`, or
  `"bic"`.

- seed:

  Optional seed to control order of term tests

## Value

For `ertte_scm_forward()` and `ertte_scm_backward()`, the updated ertte
model is returned, with the SCM history log updated internally. For
`ertte_scm_history()`, a data frame is returned containing the SCM
history log

## Details

Terms are compared with a likelihood-ratio Chi-squared test
([`stats::anova()`](https://rdrr.io/r/stats/anova.html) on nested
`survreg`/`coxph` fits) – unlike the companion `erglm` package's SCM,
there's no family-dependent choice of test here, since a
`survreg`/`coxph` model's likelihood ratio test doesn't vary by
distribution.

Three model selection criteria are available via the `criterion`
argument, mirroring the companion `emaxnls` package's development
version:

- `"p-value"` (default): a term is added if its likelihood-ratio p-value
  falls below `threshold` (forward) or removed if its p-value exceeds
  `threshold` (backward). When multiple candidates satisfy the threshold
  within a step, the one with the most extreme p-value is chosen.

- `"aic"`: a term is added (forward) or removed (backward) if doing so
  strictly decreases AIC relative to the current model. When multiple
  candidates improve AIC, the one yielding the lowest AIC is chosen.

- `"bic"`: same as `"aic"`, but using BIC as the criterion.

When `criterion` is `"aic"` or `"bic"`, the `threshold` argument has no
effect and is ignored, and `term_p_value` is left `NA` in the history
for every candidate tested that step (the likelihood-ratio test isn't
computed, since it plays no role in selection). `model_aic`/`model_bic`
are always recorded regardless of which criterion drove selection, and
the history's `criterion` column records which one was used for each
forward/backward step (`NA` for the base-model/pre-existing rows).

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

If a candidate term is aliased (perfectly collinear) with a term already
in the model, [`stats::anova()`](https://rdrr.io/r/stats/anova.html)
reports an `NA` p-value for it. That candidate is skipped for the step
(with a warning) rather than being selected or crashing the search. This
check only applies under `criterion = "p-value"`, since it's the only
criterion that computes a likelihood-ratio p-value at all – an aliased
candidate under `criterion = "aic"`/`"bic"` is simply judged (and, in
degenerate cases, potentially selected) on AIC/BIC like any other
candidate.

Two further failure modes are also handled per-candidate, rather than
aborting the whole search: if refitting with a candidate added/removed
throws an error (e.g. a single-level factor candidate that `coxph()`/
`survreg()` can't build contrasts for), that candidate is skipped with a
warning quoting the underlying error, and the rest of the candidate set
is still tried. Separately, if a candidate can't be added/removed at all
(e.g. it references a variable not present in the fitting data, so
[`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)/[`ertte_remove_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)
return `mod` unchanged), that's detected directly (the refit formula is
identical to the current model's) and the candidate is skipped with a
warning explaining why – rather than comparing the unchanged model to
itself via [`anova()`](https://rdrr.io/r/stats/anova.html), which would
produce an `NA` p-value and be misreported as aliasing/collinearity.

`candidates` is validated up front: every element must be parseable as a
formula and name exactly one covariate term (e.g. `"sex"`, not
`"sex + dose"` or `"not a formula"`).

## Examples

``` r
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
