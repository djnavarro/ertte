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

  An ertte model object, as returned by
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) or
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)

- term:

  A one-sided formula naming the term to add/remove, e.g. `~ sex`

- quiet:

  If `TRUE`, suppress the warning issued when the term can't be
  added/removed (because it's already in the model / isn't in the model,
  respectively)

## Value

An ertte model object. If the term can't be added/removed (see `quiet`),
the original `mod` is returned unchanged.

## Details

These functions are not typically called directly; they underpin
[`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
and
[`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md).
Named and shaped to match the companion `erglm` package's
`erglm_add_term()`/`erglm_remove_term()`: `term`/`candidates` are plain
formula terms, added/removed additively – categorical covariates enter
as factor levels, continuous covariates enter linearly by default or,
for a power-function parameterisation (`theta` such that
`T ~ (x / ref)^theta` on the AFT time scale, or
`h(t|x) ~ h0(t) * (x / ref)^theta` on the Cox hazard scale), by wrapping
the covariate in
[`ertte_power()`](https://ertte.djnavarro.net/reference/ertte_power.md),
e.g. `~ ertte_power(age)` or `candidates = "ertte_power(age)"`. Term
handling here works generically on formula term-labels, so
[`ertte_power()`](https://ertte.djnavarro.net/reference/ertte_power.md)
terms need no special-casing.

`mod` is refit via an internal `.ertte_refit()` helper that dispatches
on `mod`'s engine (`ertte_aft`/`ertte_coxph`) and calls the matching
constructor – so these functions (and the SCM family built on them) work
for both `ertte_aft` and `ertte_coxph` models.

## Examples

``` r
mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod2 <- ertte_add_term(mod, ~ sex)
mod3 <- ertte_remove_term(mod2, ~ sex)
```
