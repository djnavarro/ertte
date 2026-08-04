# Stepwise covariate modelling for exposure-response TTE models

Stepwise covariate modelling for exposure-response TTE models

## Usage

``` r
ertte_scm_forward(mod, candidates, threshold = 0.01, seed = NULL)

ertte_scm_backward(mod, candidates, threshold = 0.001, seed = NULL)

ertte_scm_history(mod)
```

## Arguments

- mod:

  An ertte model object

- candidates:

  Character vector with list of candidate terms

- threshold:

  Threshold to test against

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
`survreg` fits) – unlike the companion `erglm` package's SCM, there's no
family-dependent choice of test here, since a `survreg` model's
likelihood ratio test doesn't vary by distribution.

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
(with a warning) rather than being selected or crashing the search.

`candidates` is validated up front: every element must be parseable as a
formula and name exactly one covariate term (e.g. `"sex"`, not
`"sex + dose"` or `"not a formula"`).

## Examples

``` r
mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"))
ertte_scm_history(mod1)
#> # A tibble: 4 × 11
#>   iteration attempt step       action term_tested model_tested   model_converged
#>       <int>   <int> <chr>      <chr>  <chr>       <chr>          <lgl>          
#> 1         0       0 base model NA     NA          survival::Sur… TRUE           
#> 2         1       1 forward    add    ~sex        survival::Sur… TRUE           
#> 3         1       2 forward    add    ~dose       survival::Sur… TRUE           
#> 4         2       3 forward    add    ~dose       survival::Sur… TRUE           
#> # ℹ 4 more variables: term_p_value <dbl>, model_aic <dbl>, model_bic <dbl>,
#> #   model_updated <int>

mod2 <- ertte_aft(survival::Surv(time, event) ~ aucss + sex + dose, ertte_data)
mod3 <- ertte_scm_backward(mod2, candidates = c("sex", "dose"))
ertte_scm_history(mod3)
#> # A tibble: 4 × 11
#>   iteration attempt step       action term_tested model_tested   model_converged
#>       <int>   <int> <chr>      <chr>  <chr>       <chr>          <lgl>          
#> 1         0       0 base model NA     NA          survival::Sur… TRUE           
#> 2         1       1 backward   remove ~dose       survival::Sur… TRUE           
#> 3         1       2 backward   remove ~sex        survival::Sur… TRUE           
#> 4         2       3 backward   remove ~sex        survival::Sur… TRUE           
#> # ℹ 4 more variables: term_p_value <dbl>, model_aic <dbl>, model_bic <dbl>,
#> #   model_updated <int>
```
