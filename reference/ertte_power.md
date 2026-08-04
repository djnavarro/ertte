# Power-function covariate transform for exposure-response TTE models

A
[`poly()`](https://rdrr.io/r/stats/poly.html)/[`splines::ns()`](https://rdrr.io/r/splines/ns.html)-style
formula helper for entering a continuous covariate as a power function,
e.g. `ertte_power(age)` in place of a plain `age` term.

## Usage

``` r
ertte_power(x, ref = NULL)

# S3 method for class 'ertte_power'
makepredictcall(var, call)
```

## Arguments

- x:

  A numeric covariate. Every non-missing value must be strictly positive
  (`log(x / ref)` is undefined otherwise).

- ref:

  A single strictly positive reference value. Defaults to
  `median(x, na.rm = TRUE)` – the usual pop-PK/NONMEM convention of
  referencing the covariate's typical (median) value in the fitting
  data.

- var:

  The evaluated variable (here, the `ertte_power()`-transformed vector
  from the original model fit).

- call:

  The unevaluated call to be reconstructed for new data.

## Value

A numeric vector equal to `log(x / ref)`, classed `"ertte_power"`, with
`ref` stored as an attribute.

## Details

Both [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
(a log-location-scale AFT model: `log(T) = mu + scale * W`) and
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
(a Cox PH model: \`h(t \| x) = h0(t)

- exp(x'beta)`) are linear in their covariates on the model's natural (log-time or log-hazard-ratio) scale. A power-function covariate effect -- `T
  = T_ref \* (x / ref)^theta`on the AFT time scale, or`h(t \| x) = h0(t)
  \* (x /
  ref)^theta`on the Cox hazard scale -- is, after taking logs, exactly a linear term in`log(x
  / ref)\`:

    log(T)    = ... + theta * log(x / ref) + scale * W
    log(h/h0) = ... + theta * log(x / ref)

So `ertte_power(x)` reduces the power-function parameterisation to an
ordinary covariate column: the fitted `survreg()`/`coxph()` coefficient
on `ertte_power(x)` *is* the power exponent `theta` directly, and its
ordinary Wald confidence interval (from
[`confint()`](https://rdrr.io/r/stats/confint.html)/[`summary()`](https://rdrr.io/r/base/summary.html))
is exactly the confidence interval on `theta` – no delta method or
profile-likelihood machinery is needed, unlike covariate power functions
on genuinely nonlinear structural parameters (e.g. the companion
`emaxnls` package's Emax/EC50 parameters, where this reduction doesn't
apply).

`ref` is fixed at fitting time from the data `ertte_power()` is
evaluated on, and reused (not recomputed) when the fitted model is used
to predict on new data – via a `makepredictcall.ertte_power()` method,
the same mechanism
[`stats::poly()`](https://rdrr.io/r/stats/poly.html)/[`splines::ns()`](https://rdrr.io/r/splines/ns.html)
use for this purpose.

`ertte_power()` requires every non-missing value of `x` to be strictly
positive, which rules it out for covariates with a placebo/zero-dose
group (e.g. `dose`, `aucss`, `cmaxss` in `ertte_data`). This is by
design: the power-function parameterisation described in the package's
design issue is aimed at the *covariate model* (e.g. age, weight), not
the primary exposure metric, which enters the model directly (see
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)).

Nothing in
[`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)/[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/
[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
prevents combining a plain linear term (`age`) and a power term
(`ertte_power(age)`) for the same underlying variable – term handling
throughout ertte works on formula term-labels, not variable semantics,
so this is left to the user's judgement.

`makepredictcall.ertte_power()` is a
[`stats::makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
method, not typically called directly. It ensures that when a fitted
model containing an `ertte_power()` term is used to predict/simulate on
new data (via
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)/[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
on the model's [`terms()`](https://rdrr.io/r/stats/terms.html)), the
*original* fitting-time `ref` is reused rather than a new one recomputed
from whatever data is supplied – the same mechanism
[`stats::poly()`](https://rdrr.io/r/stats/poly.html)/[`splines::ns()`](https://rdrr.io/r/splines/ns.html)
use.

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
