# Landmark event probabilities in ertte

``` r

library(ertte)
library(survival)
```

This article covers
[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md),
`ertte`‘s simplest scalar exposure-response (E-R) reduction of a
time-to-event (TTE) endpoint. It’s the shorter sibling of the [RMST
article](https://ertte.djnavarro.net/articles/rmst.md), which covers a
more involved reduction (restricted mean survival time) built on the
same idea – turning a full survival curve into a single number that can
be plotted with `erplots`’ ordinary E-R grammar. If you haven’t already,
start with the [overview
article](https://ertte.djnavarro.net/articles/overview.md) for
background on
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)/
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md).

## Why a landmark?

A fitted TTE model describes an entire survival curve $`S(t \mid x)`$,
but many exposure-response questions are really about one particular
point on it: *what fraction of subjects with a given exposure profile
are expected to have had the event by some fixed, clinically meaningful
time* $`t^*`$? (“Landmark” here refers to that fixed time, not to
landmark-based dynamic prediction from repeated measurements, a
different technique with the same name.)

Reducing to a single time point converts a TTE endpoint into an ordinary
binary response, $`P(\text{event by } t^*) = 1 - S(t^* \mid x)`$, which
`erplots`’ existing `er_plot()`/`er_vpc()` grammars – built for scalar
exposure-response summaries – can visualise with no TTE-specific
plotting code at all.

## How it’s computed

[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)
doesn’t implement any new survival-curve machinery of its own. It calls
\[ertte_predict()\] at `time = landmark_time` and transforms the result:

``` math
P(\text{event by } t^*) = 1 - S(t^*).
```

Because this is a decreasing, monotonic transform, the confidence
interval bounds simply swap – the upper bound on survival becomes the
lower bound on event probability, and vice versa – with no separate
calculation needed. Whatever validity
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)’s
interval has for a given engine (a Wald interval on the linear predictor
for [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md),
or
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)’s
own log-transform interval for
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md))
carries through unchanged. This also means
[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)
inherits all of
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)’s
existing edge-case handling for free – e.g. the all-censored-Cox-model
guard, or `NA` propagation for an aliased covariate.

One consequence of delegating this way: `landmark_time` must be a single
fixed value, unlike
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)’s
`time` argument (which accepts a vector). A landmark is by definition
one specific point in time; if you want event probabilities at several
different landmarks, call
[`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)
once per landmark.

## A worked example

``` r

mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)

nd <- ertte_data[c(1, 50, 100), ]
nd[, c("id", "aucss")]
#>      id    aucss
#> 1     1 1114.089
#> 50   50  148.240
#> 100 100  427.802
```

``` r

ertte_landmark(mod_aft, nd, landmark_time = 90)
#> # A tibble: 3 × 14
#>      id sex      age weight  dose treatment aucss cmaxss event admin_censor
#>   <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl>        <dbl>
#> 1     1 Female    27     70   200 Drug      1114.  187.      1          180
#> 2    50 Male      28     58   100 Drug       148.   10.4     0          180
#> 3   100 Female    26     48   200 Drug       428.   35.1     0          180
#> # ℹ 4 more variables: landmark_time <dbl>, fit_resp <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>
```

``` r

ertte_landmark(mod_cox, nd, landmark_time = 90)
#> # A tibble: 3 × 14
#>      id sex      age weight  dose treatment aucss cmaxss event admin_censor
#>   <int> <fct>  <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl>        <dbl>
#> 1     1 Female    27     70   200 Drug      1114.  187.      1          180
#> 2    50 Male      28     58   100 Drug       148.   10.4     0          180
#> 3   100 Female    26     48   200 Drug       428.   35.1     0          180
#> # ℹ 4 more variables: landmark_time <dbl>, fit_resp <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>
```

As with
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md),
the two engines’ point estimates and intervals should generally agree
closely when both fit the data reasonably well (as here), but aren’t
expected to match exactly, since the underlying survival curves and
confidence interval constructions differ.

## Using a landmark in a visual predictive check

Passing `landmark_time` through `simulate_args` to
[`erplots::er_vpc_add_simulated()`](https://erplots.djnavarro.net/reference/er_vpc_add_simulated.html)
(or directly to `er_simulate()`/
[`simulate.ertte_model()`](https://ertte.djnavarro.net/reference/simulate.ertte_model.md))
turns each simulated replicate’s raw `sim_time`/`sim_event` into a
landmark-binary outcome for comparison against the observed data:

``` r

sim <- simulate(mod_aft, nd, nsim = 500, seed = 3157)
sim[1:3, c("id", "sim_id", "sim_time", "sim_event")]
#> # A tibble: 3 × 4
#>      id sim_id sim_time sim_event
#>   <int>  <int>    <dbl>     <dbl>
#> 1     1      1    127.          1
#> 2    50      1      0.3         0
#> 3   100      1    115.          1
```

Internally, this reduces each replicate to `1` (an event on or before
`landmark_time`), `0` (known to still be event-free at `landmark_time`),
or `NA` if the replicate was censored strictly before `landmark_time` –
a genuinely ambiguous outcome, since it’s not known whether that
replicate would have had the event by $`t^*`$ had follow-up continued.
This is the same complete-case convention used for RMST-based VPCs (see
the [RMST article](https://ertte.djnavarro.net/articles/rmst.md)’s
“Using RMST in a visual predictive check” section for the full
reasoning);
[`erplots::er_vpc_add_simulated()`](https://erplots.djnavarro.net/reference/er_vpc_add_simulated.html)’s
`mean(..., na.rm = TRUE)` aggregation excludes these `NA`s
automatically.

## Limitations

- **A landmark answers a question about one time point only.** It’s
  simpler and more directly interpretable than RMST, but it discards
  everything the fitted model says about the survival curve away from
  $`t^*`$ – two exposure profiles with very different survival curves
  can have the same landmark probability if the curves happen to cross
  near $`t^*`$.
- **The choice of $`t^*`$ matters, and isn’t estimated from the data.**
  Different landmark times can tell different (both valid) stories about
  the same fitted model, the same caveat RMST’s choice of `tau` carries.
- **A replicate censored before $`t^*`$ contributes no information to a
  landmark VPC and is dropped, not imputed or reweighted** – see “Using
  a landmark in a visual predictive check” above.

## Further reading

- The [RMST article](https://ertte.djnavarro.net/articles/rmst.md)
  covers a closely related, more involved scalar reduction (restricted
  mean survival time), including a discussion of why it needs a
  genuinely per-engine calculation, unlike
  [`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md).
- [`?ertte_landmark`](https://ertte.djnavarro.net/reference/ertte_landmark.md)
  and
  [`?ertte_predict`](https://ertte.djnavarro.net/reference/ertte_predict.md)
  for the full argument/return documentation.
