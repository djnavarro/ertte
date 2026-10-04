# Objects exported from other packages

These objects are imported from other packages. Follow the link below to
see its documentation.

## Usage

``` r
Surv(
  time,
  time2,
  event,
  type = c("right", "left", "interval", "counting", "interval2"),
  origin = 0
)
```

## Value

[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) returns an object
of class `"Surv"`, a matrix of time/event (or time/time2/event, for
interval- or counting-process-style data) values used as the response in
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)/[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
formulas. See
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) for
full details on its return value.

## Details

- survival:

  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
