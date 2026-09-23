# AFT regression modelling

Fits a parametric accelerated failure time (AFT) regression model for
time-to-event data.

## Usage

``` r
ertte_aft(formula, data, dist = "weibull", ...)
```

## Arguments

- formula:

  Model formula specifying the regression model, e.g.
  `Surv(time, event) ~ exposure`.

- data:

  Data frame containing the variables of interest.

- dist:

  The AFT distribution type to use (default is `"weibull"`).

- ...:

  Other arguments passed to
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).

## Value

A `survreg` object with additional `ertte_aft` and `ertte_model` classes
used to supply additional methods.

## Details

The `ertte_aft()` function is a thin wrapper around
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html),
used to build parametric AFT models for time-to-event data. Unlike the
original, it caches the input data frame within the returned model
object, ensuring that the data set remains accessible to downstream
tools that have access to the model object but not necessarily the
original data frame.

Like the `survreg()` function upon which it is based, `ertte_aft()` fits
a log-location-scale AFT model,

\$\$\log(t) = \mu + \sigma w\$\$

where \\\mu\\ is the linear predictor, \\\sigma\\ is the scale
parameter, and \\w\\ follows a fixed base distribution determined by
`dist`:

- `"weibull"`/`"exponential"`: \\w\\ follows a standard extreme-value
  distribution.

- `"lognormal"`: \\w\\ follows a standard normal distribution.

- `"loglogistic"`: \\w\\ follows a standard logistic distribution.

In principle `ertte_aft()` could support the full range of distributions
available to `survreg()`, but it is important to note that the extended
functionality provided by the ertte package, and hooks into other
exposure-response tools supplied by other packages such as erplots, have
not been thoroughly tested beyond four core cases. Tested values for the
`dist` argument are `weibull`, `exponential`, `lognormal`, and
`loglogistic`. It is for this reason that the output object from
`ertte_aft()` contains additional classes: the ertte-specific tools are
not necessarily valid for all `survreg` objects, only those for which it
has been explicitly defined.

Nevertheless, because the return value is genuinely a `survreg` object,
albeit with some additional information and metadata stored internally,
all the usual methods for survival regression models work unchanged,
without needing any ertte-specific equivalent. This includes
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
# fit a Weibull AFT model; because the ertte package reexports Surv(), 
# there is no need to load the survival package to build this model 
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
