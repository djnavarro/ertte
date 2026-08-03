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
  [`ertte_model()`](https://ertte.djnavarro.net/reference/ertte_model.md)

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
`erglm_add_term()`/`erglm_remove_term()`: covariates enter as plain
additive terms on the AFT location scale (linear for continuous
covariates, factor levels for categorical ones) – the richer "continuous
covariate as power function" parameterisation described in the package's
design issue is not yet implemented (see AGENTS.md).

## Examples

``` r
mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
mod2 <- ertte_add_term(mod, ~ sex)
mod3 <- ertte_remove_term(mod2, ~ sex)
```
