# Power-function covariate transform for exposure-response TTE models

A formula helper used to enter a continuous covariate as a power
transformed covariate rather than a plain linear term.

## Usage

``` r
ertte_power(x, ref = NULL)
```

## Arguments

- x:

  A numeric covariate. Every non-missing value must be strictly
  positive.

- ref:

  A single strictly positive reference value. Defaults to the median
  value of `x`.

## Value

A numeric vector equal to `log(x / ref)`, classed as an `"ertte_power"`
object, with the value of `ref` stored as an attribute.

## Details

Both [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
and
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
are linear in their covariates on the model's natural scale (i.e.,
log-time for AFT models and log-hazard-ratio for Cox models). A
power-function covariate effect is defined as follows for an AFT model:

\$\$t = t\_{\mbox{ref}} \times (x / x\_{\mbox{ref}})^\theta\$\$

The equivalent for Cox models is as follows:

\$\$h(t \| x) = h_0(t) \times (x / x\_{\mbox{ref}})^\theta\$\$

Once log-transformed, these become linear terms in `log(x / ref)`.

So `ertte_power(x)` reduces the power-function parameterisation to an
ordinary covariate column: the fitted `survreg()`/`coxph()` coefficient
on `ertte_power(x)` *is* the power exponent `theta` directly, and its
ordinary Wald confidence interval (from
[`confint()`](https://rdrr.io/r/stats/confint.html)/[`summary()`](https://rdrr.io/r/base/summary.html))
is exactly the confidence interval on `theta`.

`ertte_power()` requires every non-missing value of `x` to be strictly
positive, so cannot be used for covariates with a placebo/zero-dose
group (e.g. `dose`, `aucss`, `cmaxss` in `ertte_data`). It is not
generally appropriate for the exposure itself, only for (some) covariate
terms in the TTE model.

Nothing prevents combining a plain linear term (`age`) and a power term
(`ertte_power(age)`) for the same underlying variable in the same model
or SCM candidate set – term handling works on formula term-labels, not
variable semantics, so this is left to the user's judgement.

`ref` is fixed at fitting time from the data `ertte_power()` is
evaluated on, and reused (not recomputed) when the fitted model is used
to predict or simulate on new data, via an internal
[`stats::makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
method – the same mechanism
[`stats::poly()`](https://rdrr.io/r/stats/poly.html) uses for this
purpose.

## Examples

``` r
mod <- ertte_aft(Surv(time, event) ~ aucss + ertte_power(age), ertte_data)
summary(mod)
#> 
#> Call:
#> survival::survreg(formula = formula, data = data, dist = dist)
#>                      Value Std. Error      z       p
#> (Intercept)       4.84e+00   6.12e-02  79.19 < 2e-16
#> aucss            -6.39e-04   4.45e-05 -14.35 < 2e-16
#> ertte_power(age) -3.34e-01   2.43e-01  -1.38    0.17
#> Log(scale)       -3.37e-01   5.33e-02  -6.32 2.7e-10
#> 
#> Scale= 0.714 
#> 
#> Weibull distribution
#> Loglik(model)= -1206.3   Loglik(intercept only)= -1263.9
#>  Chisq= 115.16 on 2 degrees of freedom, p= 9.8e-26 
#> Number of Newton-Raphson Iterations: 6 
#> n= 300 
#> 

# reference value used for the power transform
attr(ertte_power(ertte_data$age), "ref")
#> [1] 27
```
