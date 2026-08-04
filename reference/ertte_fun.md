# Prediction function for an exposure-response TTE model

Prediction function for an exposure-response TTE model

## Usage

``` r
ertte_fun(object, ...)

# S3 method for class 'ertte_aft'
ertte_fun(object, ...)

# S3 method for class 'ertte_coxph'
ertte_fun(object, ...)
```

## Arguments

- object:

  An ertte model, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  (not yet implemented)

- ...:

  Passed to methods

## Value

A function with arguments `param`, `data`, and `time`.

- The `param` argument should be a vector of location coefficients;
  defaults to `coef(object)` (the fitted coefficients) if not supplied.

- The `data` argument should be a data frame or tibble; defaults to
  `object$data` (the data the model was fitted to) if not supplied.

- The `time` argument gives the time(s) at which to evaluate the
  survival function; recycled against `data`.

## Details

`ertte_fun()` is a generic, with methods for each supported engine – see
`ertte_fun.ertte_aft()`. Named `ertte_fun()` for consistency with the
companion `erglm`/`emaxnls` packages' `erglm_fun()`/`emax_fun()`, which
serve the same purpose for their respective model classes.

The `ertte_aft` method takes a fitted AFT model as input and returns a
function that evaluates the survival function `S(t)` at user-specified
parameters, data, and times (e.g. for VPCs or other counterfactual
simulation scenarios). The returned function checks that `param` is
numeric and has one entry per column of the model matrix implied by
`data`, erroring informatively rather than failing with a cryptic
"non-conformable arguments" error from matrix multiplication. `scale` is
always taken from the fitted `object`, not from `param` (the coefficient
vector from [`coef()`](https://rdrr.io/r/stats/coef.html) never includes
it).

The `ertte_coxph` method returns a function that evaluates
`S(t | x) = S0(t)^exp((x - xbar)'param)`, where `S0(t)` is the fitted
baseline survival curve (via
`survival::basehaz(object, centered = TRUE)`, held constant beyond the
last observed time, matching
[`ertte_predict.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_predict.md))
and `xbar` is `object$means` (the covariate means `coxph()` centers the
partial likelihood on when fitting – centering matters here because
`basehaz()`'s baseline is defined relative to it, not to `x = 0`). As
with `ertte_fun.ertte_aft()`, `param` only varies the linear predictor:
the baseline hazard is always taken from the fitted `object`, never
recomputed for a hypothetical `param` (that would need refitting the
partial likelihood's risk sets) – matching the level of approximation
used elsewhere in this package (e.g. `scale` for AFT models is likewise
held fixed). Since Cox models have no intercept, `param` has one entry
per covariate with no `"(Intercept)"` column, unlike
`ertte_fun.ertte_aft()`.

## Examples

``` r
mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
mod_fun <- ertte_fun(mod)

# no arguments: reproduces the fitted model's own survival predictions
s1 <- mod_fun(time = 60)

# user modifies the parameters
par2 <- coef(mod)
par2["(Intercept)"] <- par2["(Intercept)"] + 1
s2 <- mod_fun(param = par2, time = 60)

mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
mod_cox_fun <- ertte_fun(mod_cox)

# no arguments: reproduces the fitted model's own survival predictions
s1 <- mod_cox_fun(time = 60)

# user modifies the parameters
par2 <- coef(mod_cox)
par2["aucss"] <- par2["aucss"] * 1.5
s2 <- mod_cox_fun(param = par2, time = 60)
```
