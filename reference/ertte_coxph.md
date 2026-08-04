# Fit an exposure-response time-to-event Cox PH model based on `coxph()`

Fit an exposure-response time-to-event Cox PH model based on `coxph()`

## Usage

``` r
ertte_coxph(formula, data, ...)
```

## Arguments

- formula:

  Model formula, with a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  object as the response, e.g. `Surv(time, event) ~ exposure`.

- data:

  Data set

- ...:

  Other arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## Value

A coxph object with extra `ertte_coxph`/`ertte_model` classes

## Details

The returned object has class
`c("ertte_coxph", "ertte_model", "coxph")`: it *is* a `coxph` object,
with a little extra metadata attached. This means all of the usual
`coxph` methods work unchanged, without needing an ertte-specific
equivalent – e.g. [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html), and
[`anova()`](https://rdrr.io/r/stats/anova.html) for comparing nested
models.

`ertte_coxph()` is the semi-parametric sibling of
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md).
Unlike
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md),
there's no `dist` argument: Cox PH doesn't assume a parametric baseline
hazard, so there's no distribution to select. Both share the
`"ertte_model"` superclass, so functions that only need generic
operations ([`update()`](https://rdrr.io/r/stats/update.html),
[`anova()`](https://rdrr.io/r/stats/anova.html), the SCM family) work
across either engine; functions with engine-specific logic (e.g.
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md),
[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md))
dispatch via the `"ertte_aft"`/`"ertte_coxph"` subclass.

[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
and [`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md)
have `ertte_coxph` methods (see
[`ertte_predict.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_predict.md)/[`ertte_fun.ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_fun.md)),
both built on the fitted baseline hazard.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) works too, via the
shared
[`simulate.ertte_model()`](https://ertte.djnavarro.net/reference/simulate.ertte_model.md)
method – no separate `simulate.ertte_coxph()` is needed, since it
dispatches internally (via `.ertte_simulate_draws()`) on engine.

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
