# Add or remove a covariate term from an exposure-response TTE model

Add or remove a single covariate term from an existing ertte model,
returning a new fitted model object.

## Usage

``` r
ertte_add_term(mod, term, quiet = FALSE)

ertte_remove_term(mod, term, quiet = FALSE)
```

## Arguments

- mod:

  An ertte model object.

- term:

  A one-sided formula naming the term to add/remove, e.g. `~ sex`.

- quiet:

  Should warnings be suppressed? Defaults to `FALSE`.

## Value

An ertte model object. If the term can't be added/removed (see `quiet`),
the original `mod` is returned unchanged.

## Details

These functions are not typically called directly; they underpin
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
and
[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md),
used to add or remove a single term from a time-to-event regression
model. Regardless of whether the model is a parametric AFT model or a
Cox proportional hazard model, the `term` to be added or removed is
defined by a single one-sided formula. Categorical covariates enter the
model as factor levels, whereas continuous covariates enter as linear
terms by default. However, for a power-function parameterisation
(`theta` such that `T ~ (x /ref)^theta` on the AFT time scale, or
`h(t|x) ~ h0(t) * (x / ref)^theta` on the Cox hazard scale), the
covariate can be wrapped in
[`ertte_power()`](https://ertte.djnavarro.net/reference/ertte_power.md),
e.g. `~ ertte_power(age)`.

## Examples

``` r
mod1 <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod2 <- ertte_add_term(mod1, ~ sex)
mod3 <- ertte_remove_term(mod2, ~ sex)
```
