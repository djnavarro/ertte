# Stepwise covariate modelling in ertte

``` r

library(ertte)
library(survival)
```

This article covers `ertte`’s stepwise covariate modelling (SCM)
machinery:
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md),
the single-term
[`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)/[`ertte_remove_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)
helpers they’re built on, and
[`ertte_scm_history()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
for inspecting what happened during a search. It assumes familiarity
with
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)/[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
– see the [overview
article](https://ertte.djnavarro.net/articles/overview.md) if you
haven’t already read it.

## What SCM does

Stepwise covariate modelling is a common approach (particularly in
pharmacometrics) for deciding which covariates belong in a model:
starting from a base model, candidate covariate terms are added (forward
selection) or removed (backward elimination) one at a time, keeping a
term only if it improves the model fit by more than some threshold
amount.

`ertte` implements this as
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md),
working identically for
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) and
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
models – both functions dispatch internally on the model’s engine, so no
separate AFT/Cox versions are needed at the user-facing level.

## Forward selection

`ertte_scm_forward(mod, candidates, threshold)` starts from a fitted
model and repeatedly tests each not-yet-included term in `candidates`,
adding whichever single term improves the fit the most (by
likelihood-ratio test p-value) at each step, provided that improvement
clears `threshold`. It stops when no remaining candidate clears the
threshold.

``` r

mod0 <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod1 <- ertte_scm_forward(
  mod0,
  candidates = c("sex", "dose", "age", "weight"),
  threshold = 0.01,
  seed = 4821
)
coef(mod1)
#>   (Intercept)         aucss       sexMale 
#>  4.7360379087 -0.0006680892  0.2743522150
```

Every step (including unsuccessful ones) is logged, and can be inspected
with
[`ertte_scm_history()`](https://ertte.djnavarro.net/reference/ertte_scm.md):

``` r

ertte_scm_history(mod1)
#> # A tibble: 8 × 12
#>   iteration attempt step       criterion action term_tested model_tested        
#>       <int>   <int> <chr>      <chr>     <chr>  <chr>       <chr>               
#> 1         0       0 base model NA        NA     NA          Surv(time, event) ~…
#> 2         1       1 forward    p-value   add    ~sex        Surv(time, event) ~…
#> 3         1       2 forward    p-value   add    ~weight     Surv(time, event) ~…
#> 4         1       3 forward    p-value   add    ~dose       Surv(time, event) ~…
#> 5         1       4 forward    p-value   add    ~age        Surv(time, event) ~…
#> 6         2       6 forward    p-value   add    ~dose       Surv(time, event) ~…
#> 7         2       7 forward    p-value   add    ~weight     Surv(time, event) ~…
#> 8         2       8 forward    p-value   add    ~age        Surv(time, event) ~…
#> # ℹ 5 more variables: model_converged <lgl>, term_p_value <dbl>,
#> #   model_aic <dbl>, model_bic <dbl>, model_updated <int>
```

Each row is one *attempt* – one candidate term tested within one
iteration of the search. `term_p_value` is the likelihood-ratio test
p-value comparing the model with and without that term; `model_updated`
is `1` for the single attempt (if any) that was actually kept at the end
of that iteration, and `0` for the rest. An iteration with every
`model_updated == 0` means no candidate cleared the threshold that
round, which is also when the search stops.

## Backward elimination

`ertte_scm_backward(mod, candidates, threshold)` works the other
direction: starting from a model that already contains the candidate
terms, it repeatedly tests removing each one, dropping whichever single
term’s removal is *least* damaging (highest likelihood-ratio p-value),
provided that p-value clears `threshold`. Because removing a term should
be a more conservative decision than adding one,
[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)’s
default `threshold` (`0.001`) is stricter than
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)’s
(`0.01`) – the usual convention in pharmacometric SCM workflows,
sometimes run together as “forward inclusion at p \< 0.01, backward
elimination at p \< 0.001”.

``` r

mod_full <- ertte_aft(Surv(time, event) ~ aucss + sex + dose + age, ertte_data)
mod_back <- ertte_scm_backward(
  mod_full,
  candidates = c("sex", "dose", "age"),
  threshold = 0.001,
  seed = 6039
)
coef(mod_back)
#>   (Intercept)         aucss 
#>  4.8563375087 -0.0006407913
```

## The significance test

Every comparison in `ertte`’s SCM machinery is a likelihood-ratio
Chi-squared test, via
[`stats::anova()`](https://rdrr.io/r/stats/anova.html) on two nested
`survreg`/`coxph` fits. Unlike the companion `erglm` package (whose
GLM-based SCM needs a family-dependent choice of test), there’s no
equivalent choice to make here: a `survreg`/`coxph` model’s
likelihood-ratio test doesn’t vary by distribution or engine.

If a candidate term is aliased (perfectly collinear) with a term already
in the model, [`anova()`](https://rdrr.io/r/stats/anova.html) reports an
`NA` p-value for the comparison. Rather than crash the search or
silently select/reject the term,
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
skip that candidate for the current step (with a warning).

## Selection criteria: p-value, AIC, or BIC

By default, both
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
and
[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
select terms using the likelihood-ratio p-value described above,
compared against `threshold`. The `criterion` argument lets you switch
to an information-criterion-based rule instead:

``` r

mod_aic <- ertte_scm_forward(
  mod0,
  candidates = c("sex", "dose", "age", "weight"),
  criterion = "aic",
  seed = 4821
)
coef(mod_aic)
#>   (Intercept)         aucss       sexMale 
#>  4.7360379087 -0.0006680892  0.2743522150
```

``` r

mod_bic <- ertte_scm_backward(
  mod_full,
  candidates = c("sex", "dose", "age"),
  criterion = "bic",
  seed = 6039
)
coef(mod_bic)
#>   (Intercept)         aucss       sexMale 
#>  4.7360379087 -0.0006680892  0.2743522150
```

With `criterion = "aic"` or `criterion = "bic"`, a term is added
(forward) or removed (backward) only if doing so strictly decreases the
chosen information criterion relative to the current model – if several
candidates would each improve it, the one giving the lowest resulting IC
is kept. `threshold` has no effect in this mode and is silently ignored.

The history still records everything, via a new `criterion` column that
shows which rule was in force for each step; `term_p_value` is left `NA`
for IC-driven steps, since the likelihood-ratio test isn’t computed at
all when it plays no role in selection (`model_aic`/`model_bic` are
always populated regardless of `criterion`, as they always were):

``` r

ertte_scm_history(mod_aic)
#> # A tibble: 8 × 12
#>   iteration attempt step       criterion action term_tested model_tested        
#>       <int>   <int> <chr>      <chr>     <chr>  <chr>       <chr>               
#> 1         0       0 base model NA        NA     NA          Surv(time, event) ~…
#> 2         1       1 forward    aic       add    ~sex        Surv(time, event) ~…
#> 3         1       2 forward    aic       add    ~weight     Surv(time, event) ~…
#> 4         1       3 forward    aic       add    ~dose       Surv(time, event) ~…
#> 5         1       4 forward    aic       add    ~age        Surv(time, event) ~…
#> 6         2       6 forward    aic       add    ~dose       Surv(time, event) ~…
#> 7         2       7 forward    aic       add    ~weight     Surv(time, event) ~…
#> 8         2       8 forward    aic       add    ~age        Surv(time, event) ~…
#> # ℹ 5 more variables: model_converged <lgl>, term_p_value <dbl>,
#> #   model_aic <dbl>, model_bic <dbl>, model_updated <int>
```

## Adding or removing a single term directly

[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
are both built on
[`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)/[`ertte_remove_term()`](https://ertte.djnavarro.net/reference/ertte_term.md),
which are also exported and useful directly when you want to test one
specific term without running a full search:

``` r

mod2 <- ertte_add_term(mod0, ~sex)
anova(mod0, mod2)
#>         Terms Resid. Df    -2*LL Test Df Deviance    Pr(>Chi)
#> 1       aucss       297 2414.530      NA       NA          NA
#> 2 aucss + sex       296 2406.311 +sex  1 8.219176 0.004145002
```

``` r

mod3 <- ertte_remove_term(mod2, ~sex)
identical(coef(mod0), coef(mod3))
#> [1] TRUE
```

`term` is a one-sided formula naming exactly one term – categorical
covariates enter as factor levels, continuous covariates enter linearly
by default. Refitting happens via an internal helper that dispatches on
the model’s engine and calls the matching constructor
([`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)/
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md))
directly, rather than
[`stats::update()`](https://rdrr.io/r/stats/update.html) –
[`update()`](https://rdrr.io/r/stats/update.html) doesn’t work on
`ertte` model objects, because the fitted object’s `$call` refers to the
constructor’s own local `formula`/`data` bindings, not anything visible
in the caller’s frame. This is transparent to normal usage; it only
matters if you were tempted to call
[`update()`](https://rdrr.io/r/stats/update.html) on an `ertte` model
yourself.

## Power-function covariates

Continuous covariates don’t have to enter a candidate set linearly.
`ertte_power(x)` (see
[`?ertte_power`](https://ertte.djnavarro.net/reference/ertte_power.md)
for the full derivation) reparametrises a continuous covariate as a
power function – `T = T_ref * (x / ref)^theta` on the AFT time scale, or
`h(t | x) = h0(t) * (x / ref)^theta` on the Cox hazard scale – which,
after taking logs, is exactly a linear term in `log(x / ref)`. This
means it slots into
[`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)/
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
with no special handling: term handling throughout `ertte` works
generically on formula term-labels, so `~ ertte_power(age)` (or
`"ertte_power(age)"` in a `candidates` vector) is already a legitimate
term/candidate.

``` r

mod4 <- ertte_add_term(mod0, ~ ertte_power(age))
coef(mod4)["ertte_power(age)"]
#> ertte_power(age) 
#>       -0.3342985
```

The fitted coefficient on `ertte_power(age)` *is* the power exponent
$`\theta`$ directly, and its ordinary Wald confidence interval (from
[`confint()`](https://rdrr.io/r/stats/confint.html)) is already the
confidence interval on $`\theta`$ – no delta method or profile
likelihood needed.
[`ertte_power()`](https://ertte.djnavarro.net/reference/ertte_power.md)
requires every non-missing value of its input to be strictly positive,
which is why it’s suited to covariates like age or weight rather than an
exposure metric with a placebo/zero-dose group (e.g. `aucss`, which has
`0`s for the placebo arm in `ertte_data`).

A candidate set can freely mix plain linear terms and power terms for
different (or even the same) underlying variable – `ertte` doesn’t
reason about variable semantics, so combining `age` and
`ertte_power(age)` in the same candidate set, or the same model, is left
to the analyst’s judgement.

## Summary of conventions

- **The default `criterion = "p-value"` uses thresholds on a
  likelihood-ratio Chi-squared test** – this matches the classic
  pharmacometric SCM procedure. Set `criterion = "aic"` or `"bic"` to
  select on an information criterion instead, in which case `threshold`
  is ignored.
  [`ertte_scm_history()`](https://ertte.djnavarro.net/reference/ertte_scm.md)’s
  `model_aic`/`model_bic` columns are always logged at every step
  regardless of `criterion`, so you can inspect the search through an
  information-criterion lens after the fact even when it wasn’t what
  drove selection.
- **`seed` only affects the order candidates are tested within a step**,
  which only matters in the essentially measure-zero case of an exact
  p-value tie between competing candidates. Model fitting itself is
  deterministic given a formula and dataset.
- **An aliased candidate (`NA` p-value) is skipped with a warning**, not
  treated as either a pass or a failure.
- **Both engines share the same SCM code path** via the internal
  `.ertte_refit()` dispatch described above – there’s nothing
  engine-specific to configure when switching between
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
  and
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  models.

## Further reading

- [`?ertte_scm`](https://ertte.djnavarro.net/reference/ertte_scm.md) and
  [`?ertte_term`](https://ertte.djnavarro.net/reference/ertte_term.md)
  for the full argument/return documentation.
- [`?ertte_power`](https://ertte.djnavarro.net/reference/ertte_power.md)
  for the derivation behind the power-function parameterisation.
- The [overview
  article](https://ertte.djnavarro.net/articles/overview.md) for
  background on the AFT/Cox PH model families these functions build on.
