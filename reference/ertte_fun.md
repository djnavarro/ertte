# Prediction function for an exposure-response TTE model

Returns a function that evaluates a fitted ertte model's survival
function at user-specified data, times, and (optionally) counterfactual
parameters, without needing to refit the model.

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

  An ertte model.

- ...:

  Passed to methods.

## Value

A function with arguments `data`, `time`, and `param`:

- The `data` argument should be a data frame or tibble; defaults to
  `object$data` (the data the model was fitted to) if not supplied.

- The `time` argument gives the time(s) at which to evaluate the
  survival function; recycled against `data`.

- The `param` argument should be a vector of location coefficients;
  defaults to `coef(object)` (the fitted coefficients) if not supplied.

The function returns a vector of survival probabilities.

## Details

`ertte_fun()` is a generic function, with methods for each of the
supported exposure-response time-to-event model classes.

## AFT models

The `ertte_aft` method takes a fitted AFT model as input and returns a
function that evaluates the survival function `S(t)` at user-specified
parameters, data, and times. Note that the `scale` parameter is always
taken from the fitted object, and not passed via `param`.

## Cox PH models

The `ertte_coxph` method takes a fitted Cox model as input and similarly
returns a function that evaluates the survival probabilities. More
precisely, it returns

\$\$S(t\|x) = S_0(t)^{\exp((x - \bar{x})' \beta)}\$\$

where \\x\\ is the vector of covariates, \\\bar{x}\\ is the mean of the
covariates in the fitted model, \\\beta\\ is the vector of coefficients,
and \\(x - \bar{x})' \beta\\ (a matrix multiplication) is the linear
predictor that \\S_0(t)\\ is raised to the power of.

The baseline survival function \\S_0(t)\\ refers to the fitted baseline
survival curve, calulated by calling
[`survival::basehaz()`](https://rdrr.io/pkg/survival/man/basehaz.html)
with `centered = TRUE` and held constant beyond the last observed time,
consistent with the approach adopted by
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
for Cox proportional hazard models.

Analogous to the AFT model case, `param` is only used to vary the linear
predictor: the baseline hazard is always taken from the fitted model,
and not recomputed for a hypothetical `param` value as that would
require refitting the model.

## Examples

``` r
# fit an AFT model
mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod_aft_fun <- ertte_fun(mod_aft)

# when called with no arguments it reproduces the fitted 
# model's own survival predictions
s1 <- mod_aft_fun(time = 60)
s1[1:5]
#> [1] 0.3925581 0.5653608 0.7080801 0.7080801 0.2939277

# when called with an explicit parameter vector, the model
# is evaluated at those modified parameters
par_new <- coef(mod_aft)
par_new["(Intercept)"] <- par_new["(Intercept)"] + 1
s2 <- mod_aft_fun(param = par_new, time = 60)
s2[1:5]
#> [1] 0.7933096 0.8683024 0.9180736 0.9180736 0.7384583

# the same logic applies for Cox regression models 
mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
mod_cox_fun <- ertte_fun(mod_cox)
s3 <- mod_cox_fun(time = 60)
s3[1:5]
#> [1] 0.3752982 0.5485030 0.6940022 0.6940022 0.2780437
```
