# Fit an exposure-response time-to-event model based on `survreg()`

Fit an exposure-response time-to-event model based on `survreg()`

## Usage

``` r
ertte_model(formula, data, dist = "weibull", ...)
```

## Arguments

- formula:

  Model formula, with a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  object as the response, e.g. `Surv(time, event) ~ exposure`.

- data:

  Data set

- dist:

  The AFT distribution to fit, as for
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).
  Defaults to `"weibull"`. Tested and officially supported for
  `"exponential"`, `"weibull"`, `"lognormal"`, and `"loglogistic"` – see
  [`ertte_select_distribution()`](https://ertte.djnavarro.net/reference/ertte_select_distribution.md)
  for choosing among them by AIC.

- ...:

  Other arguments passed to
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).

## Value

A survreg object with an extra `ertte_model` class

## Details

The returned object has class `c("ertte_model", "survreg")`: it *is* a
`survreg` object, with a little extra metadata attached. This means all
of the usual `survreg` methods work unchanged, without needing an
ertte-specific equivalent – e.g.
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html), and
[`anova()`](https://rdrr.io/r/stats/anova.html) for comparing nested
models.
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
is a separate, ertte-specific alternative to
[`predict()`](https://rdrr.io/r/stats/predict.html) that returns
survival probabilities with confidence intervals in a tidy data frame;
the two are complementary, not competing.

All four supported distributions are log-location-scale AFT models:
`log(T) = mu + scale * W`, where `mu` is the linear predictor
(intercept + covariates) and `W` follows a distribution that depends
only on `dist` (extreme-value for `"exponential"`/`"weibull"`, standard
normal for `"lognormal"`, standard logistic for `"loglogistic"`) – see
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
and `.ertte_dist_info()`.

## Examples

``` r
mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
mod
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

# other AFT distributions are also supported
mod_ln <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data, dist = "lognormal")
mod_ln
#> Call:
#> survival::survreg(formula = formula, data = data, dist = dist)
#> 
#> Coefficients:
#>   (Intercept)         aucss 
#>  4.4908037548 -0.0006546191 
#> 
#> Scale= 0.978441 
#> 
#> Loglik(model)= -1222.3   Loglik(intercept only)= -1271.4
#>  Chisq= 98.21 on 1 degrees of freedom, p= <2e-16 
#> n= 300 
```
