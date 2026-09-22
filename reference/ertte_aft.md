# AFT regression modeling

Fits a parametric accelerated failure time (AFT) regression model for
time-to-event data via
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).

## Usage

``` r
ertte_aft(formula, data, dist = "weibull", ...)
```

## Arguments

- formula:

  Model formula specifying the regression model, e.g.
  `Surv(time, event) ~ exposure`.

- data:

  Data set containing the variables of interest.

- dist:

  The AFT distribution type to use, defaulting to `"weibull"`.

- ...:

  Other arguments passed to
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).

## Value

A `survreg` object with additional `ertte_aft` and `ertte_model` classes
used to supply additional methods.

## Details

Like the `survreg()` function upon which it is based, `ertte_aft()`
supports four log-location-scale AFT models of the form
`log(t) = mu + scale * w`, where `mu` is the linear predictor and the
distribution of `w` is dependent on the choice of `dist`: extreme-value
distributions for `"exponential"` and `"weibull"` models, a standard
normal for `"lognormal"`, and a standard logistic for `"loglogistic"`.

Because the return value is a `survreg` object, all the usual methods
for survival regression models work unchanged, without neeing an
ertte-specific equivalent. This includes
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html), and
[`anova()`](https://rdrr.io/r/stats/anova.html). Additional methods
supplied via the ertte-specific classes include
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md),
[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md), and
[`simulate()`](https://rdrr.io/r/stats/simulate.html).

## Examples

``` r
# fit a Weibull AFT model
mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
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
mod_ln <- ertte_aft(Surv(time, event) ~ aucss, ertte_data, dist = "lognormal")
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
