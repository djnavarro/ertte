# Cox proportional hazard modeling

Fits a semi-parametric Cox proportional-hazards regression of
time-to-event on covariates via
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## Usage

``` r
ertte_coxph(formula, data, ...)
```

## Arguments

- formula:

  Model formula specifying the regression model, e.g.
  `Surv(time, event) ~ exposure`.

- data:

  Data set containing the variables of interest.

- ...:

  Other arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## Value

A coxph object with with additional `ertte_coxph` and `ertte_model`
classes used to supply additional methods.

## Details

Because the return value is a `coxph` object, all the usual methods for
survival regression models work unchanged, without neeing an
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
mod <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
mod
#> Call:
#> survival::coxph(formula = formula, data = data, model = TRUE)
#> 
#>            coef exp(coef)  se(coef)     z      p
#> aucss 8.859e-04 1.001e+00 7.515e-05 11.79 <2e-16
#> 
#> Likelihood ratio test=110.2  on 1 df, p=< 2.2e-16
#> n= 300, number of events= 232 
```
