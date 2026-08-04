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

## Not yet implemented

[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md),
[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md), and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) don't yet have
`ertte_coxph` methods – unlike the closed-form survival function
available for AFT models, Cox PH prediction/simulation needs a baseline
hazard estimate (e.g. via
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)),
which is separate follow-up work (see AGENTS.md). Calling any of these
on an `ertte_coxph` object currently errors with "no applicable method".

## Examples

``` r
mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
mod
#> Call:
#> survival::coxph(formula = formula, data = data)
#> 
#>            coef exp(coef)  se(coef)     z      p
#> aucss 8.859e-04 1.001e+00 7.515e-05 11.79 <2e-16
#> 
#> Likelihood ratio test=110.2  on 1 df, p=< 2.2e-16
#> n= 300, number of events= 232 
```
