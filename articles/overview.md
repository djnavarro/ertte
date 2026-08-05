# Survival models in ertte

``` r

library(ertte)
library(survival)
```

This article is the place to start if you’re new to `ertte`. It walks
through what a time-to-event (TTE) endpoint looks like, a short
refresher on the two statistical ideas `ertte`’s model engines are built
on (parametric accelerated failure time models and the Cox
proportional-hazards model), and how to fit, inspect, and predict from
both engines via
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) and
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md).
It assumes general familiarity with regression modelling, but not with
survival analysis specifically.

Three other articles build on this one: [restricted mean survival
time](https://ertte.djnavarro.net/articles/rmst.md) and [landmark event
probabilities](https://ertte.djnavarro.net/articles/landmark.md) cover
`ertte`’s two ways of reducing a TTE endpoint to a single scalar
exposure-response summary, and [stepwise covariate
modelling](https://ertte.djnavarro.net/articles/scm.md) covers `ertte`’s
tools for building up a covariate model.

## Time-to-event data

A time-to-event outcome records, for each subject, how long it took for
some event of interest to happen – and, crucially, whether it was
actually observed to happen at all. `ertte_data` (a simulated dataset
shipped with the package) has the shape this kind of data usually takes:

``` r

ertte_data[1:6, c("id", "aucss", "time", "event")]
#>   id    aucss  time event
#> 1  1 1114.089  77.4     1
#> 2  2  561.267  26.8     0
#> 3  3    0.000 180.0     0
#> 4  4    0.000  16.8     0
#> 5  5 1415.503  33.9     1
#> 6  6    0.000  81.3     1
```

`time` is the number of days until either the event occurred or the
subject was last observed; `event` is `1` if the event was actually
observed by `time`, and `0` if the subject was instead **censored** –
still event-free the last time they were seen, but with no way to know
what happened to them afterwards (they left the study, the study ended,
etc.). Both engines in `ertte` take a formula response built with
`survival::Surv(time, event)`, which encodes exactly this: an observed
time, plus whether it represents a true event or a censoring.

Handling censored observations correctly – using the information that a
censored subject survived at least that long, without pretending to know
what happened after – is the entire reason TTE models exist as a
distinct class of model, rather than just regressing `time` directly.

## A quick refresher: the survival and hazard functions

Two closely related functions describe a TTE outcome’s distribution.

The **survival function** $`S(t) = P(T > t)`$ is the probability of
still being event-free at time $`t`$. It starts at $`S(0) = 1`$ and
decreases (weakly) as $`t`$ increases.

The **hazard function** $`h(t)`$ is the instantaneous event rate at time
$`t`$, conditional on having survived to $`t`$:

``` math
h(t) = \lim_{\Delta t \to 0} \frac{P(t \le T < t + \Delta t \mid T \ge t)}{\Delta t}.
```

The two are related by $`S(t) = \exp\left(-\int_0^t h(u)\, du\right)`$:
the survival function is fully determined by the hazard, and vice versa.
Every model discussed below is really a model for one of these two
functions (optionally as a function of covariates $`x`$) –
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
models $`S(t
\mid x)`$ directly, and
[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
models $`h(t \mid x)`$ directly, and each can recover the other.

## Two engines, two ways of modelling $`S(t \mid x)`$

`ertte` fits two structurally different kinds of TTE model. Both accept
the same kind of formula and the same kind of covariates, but make
different assumptions and answer slightly different questions.

### `ertte_aft()`: parametric accelerated failure time models

An accelerated failure time (AFT) model assumes covariates act
multiplicatively on the *time scale* – a covariate can speed up or slow
down the passage toward the event, like a “clock” running fast or slow.
Equivalently, on the log-time scale, it’s an ordinary linear model:

``` math
\log(T) = \mu(x) + \sigma W,
```

where $`\mu(x)`$ is the linear predictor (intercept + covariate
effects), $`\sigma`$ is a scale parameter, and $`W`$ is a random
variable following a fixed “base” distribution that doesn’t depend on
$`x`$. Different choices of base distribution give different named AFT
models, all supported by
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)’s
`dist` argument:

| `dist` | Base distribution of $`W`$ | Named model |
|----|----|----|
| `"exponential"` | Standard extreme-value, fixed $`\sigma = 1`$ | Exponential |
| `"weibull"` (default) | Standard extreme-value | Weibull |
| `"lognormal"` | Standard normal | Log-normal |
| `"loglogistic"` | Standard logistic | Log-logistic |

Because the whole family shares this log-location-scale structure, all
four give a fully parametric, closed-form survival function

``` math
S(t \mid x) = 1 - F\!\left(\frac{\log t - \mu(x)}{\sigma}\right),
```

where $`F`$ is the base distribution’s CDF.
[`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
wraps
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html),
and adds two `ertte_aft`/`ertte_model` classes on top so it
interoperates with the rest of the package (and, if installed, with
`erplots`).

``` r

mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
mod_aft
#> Call:
#> survival::survreg(formula = formula, data = data, dist = dist)
#> 
#> Coefficients:
#>   (Intercept)         aucss 
#>  4.8563375087 -0.0006407913 
#> 
#> Scale= 0.7164032 
#> 
#> Loglik(model)= -1207.3   Loglik(intercept only)= -1263.9
#>  Chisq= 113.25 on 1 degrees of freedom, p= <2e-16 
#> n= 300
```

Coefficients are on the log-time scale, so `exp(coef(...))` is usually
the more interpretable quantity: a **time ratio**, i.e. the
multiplicative change in expected survival time per unit change in the
covariate.

``` r

exp(stats::coef(mod_aft))
#> (Intercept)       aucss 
#> 128.5525164   0.9993594
```

Here, `aucss`’s time ratio is just under 1: higher exposure is
associated with a (slightly) *shorter* time to event, holding nothing
else fixed – consistent with `ertte_data`’s known Weibull AFT ground
truth (see
[`?ertte_data`](https://ertte.djnavarro.net/reference/ertte_data.md)).

Choosing among the four supported distributions by AIC is automated by
[`ertte_aft_select_distribution()`](https://ertte.djnavarro.net/reference/ertte_aft_select_distribution.md):

``` r

cmp <- ertte_aft_select_distribution(Surv(time, event) ~ aucss, ertte_data)
cmp$comparison
#> # A tibble: 4 × 5
#>   dist        logLik   aic   bic converged
#>   <chr>        <dbl> <dbl> <dbl> <lgl>    
#> 1 weibull     -1207. 2421. 2432. TRUE     
#> 2 loglogistic -1216. 2438. 2449. TRUE     
#> 3 lognormal   -1222. 2451. 2462. TRUE     
#> 4 exponential -1224. 2452. 2460. TRUE
```

### `ertte_coxph()`: the semi-parametric Cox proportional-hazards model

The Cox model instead assumes covariates act multiplicatively on the
*hazard scale*, and – unlike the AFT family – makes no assumption about
the shape of the baseline hazard over time:

``` math
h(t \mid x) = h_0(t) \exp(x'\beta),
```

where $`h_0(t)`$ is an unspecified **baseline hazard** (the hazard for a
subject with $`x = 0`$, or more precisely $`x`$ at the values
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) centers on) and
$`\beta`$ are the log-hazard-ratio coefficients. This is a
*semi-parametric* model: $`\beta`$ is estimated by maximising a partial
likelihood that doesn’t require ever specifying $`h_0(t)`$’s functional
form, and $`h_0(t)`$ itself (if needed) is estimated afterwards,
nonparametrically, from the same fit (via
[`survival::basehaz()`](https://rdrr.io/pkg/survival/man/basehaz.html),
Breslow’s estimator).

The name “proportional hazards” describes its central, testable
assumption: the *ratio* of hazards between any two covariate profiles,
$`h(t \mid
x_1)/h(t \mid x_2) = \exp((x_1 - x_2)'\beta)`$, is constant over time –
it doesn’t depend on $`t`$. This is a genuinely different assumption
from the AFT family’s “constant multiplicative effect on the time
scale”, and the two coincide only for the Weibull/exponential AFT model
(the only member of the AFT family that is *also* a proportional-hazards
model).

``` r

mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
mod_cox
#> Call:
#> survival::coxph(formula = formula, data = data, model = TRUE)
#> 
#>            coef exp(coef)  se(coef)     z      p
#> aucss 8.859e-04 1.001e+00 7.515e-05 11.79 <2e-16
#> 
#> Likelihood ratio test=110.2  on 1 df, p=< 2.2e-16
#> n= 300, number of events= 232
```

Coefficients are on the log-hazard-ratio scale, so `exp(coef(...))`
gives a **hazard ratio**: the multiplicative change in the instantaneous
event rate per unit change in the covariate, at any time $`t`$.

``` r

exp(stats::coef(mod_cox))
#>    aucss 
#> 1.000886
```

A hazard ratio just above 1 here says higher `aucss` is associated with
a slightly *higher* instantaneous event rate at any given time – the
same underlying exposure effect as the AFT model above, described on a
different scale (a time ratio below 1 and a hazard ratio above 1 are two
consistent ways of saying “higher exposure, worse outcome”).

The proportional-hazards assumption itself can be checked with
[`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html),
which tests whether a covariate’s effect on the hazard actually stays
constant over time (a significant result suggests it doesn’t):

``` r

cox.zph(mod_cox)
#>        chisq df    p
#> aucss  0.574  1 0.45
#> GLOBAL 0.574  1 0.45
```

Here there’s no evidence against the assumption for `aucss` –
unsurprising given `ertte_data` was simulated from a Weibull AFT ground
truth, which (as noted above) is also a proportional-hazards model.

### Choosing between them

Both engines answer closely related exposure-response questions, and for
data that’s well described by a Weibull model, their conclusions should
agree in substance (as above). Some practical differences worth knowing
about:

- [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)
  gives a fully parametric survival curve defined for all $`t
  > 0`$, and a `dist`-dependent closed form for effects like time ratios
  – useful when extrapolating beyond the observed follow-up range, or
  when a parametric quantity like RMST needs integrating over the whole
  curve (see the [RMST
  article](https://ertte.djnavarro.net/articles/rmst.md)).
- [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  makes no assumption about the *shape* of the baseline hazard over
  time, which is often a safer default when there’s no strong reason to
  believe a particular parametric family. Its cost: no closed-form
  survival curve beyond the last observed follow-up time
  ([`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)/[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md)
  hold survival flat past that point, rather than extrapolating a
  curve), and no way to interpret effects as time ratios.

Both share the same `"ertte_model"` superclass, so most of the package’s
downstream machinery –
[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md),
[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md),
[`simulate()`](https://rdrr.io/r/stats/simulate.html), and the stepwise
covariate modelling functions – work identically regardless of which
engine a model came from; see “API naming: AFT vs Cox PH” in the
package’s development notes for the full dispatch scheme.

## Predicting and simulating

[`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
returns survival probabilities (with confidence intervals) at one or
more times, for one or more covariate profiles, for either engine:

``` r

nd <- ertte_data[c(1, 50, 100), ]
ertte_predict(mod_aft, nd, time = c(90, 180))
#> # A tibble: 6 × 14
#>      id sex     age weight  dose treatment aucss cmaxss  time event admin_censor
#>   <int> <fct> <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>        <dbl>
#> 1     1 Fema…    27     70   200 Drug      1114.  187.     90     1          180
#> 2     1 Fema…    27     70   200 Drug      1114.  187.    180     1          180
#> 3    50 Male     28     58   100 Drug       148.   10.4    90     0          180
#> 4    50 Male     28     58   100 Drug       148.   10.4   180     0          180
#> 5   100 Fema…    26     48   200 Drug       428.   35.1    90     0          180
#> 6   100 Fema…    26     48   200 Drug       428.   35.1   180     0          180
#> # ℹ 3 more variables: fit_survival <dbl>, ci_lower <dbl>, ci_upper <dbl>
ertte_predict(mod_cox, nd, time = c(90, 180))
#> # A tibble: 6 × 14
#>      id sex     age weight  dose treatment aucss cmaxss  time event admin_censor
#>   <int> <fct> <int>  <dbl> <dbl> <fct>     <dbl>  <dbl> <dbl> <dbl>        <dbl>
#> 1     1 Fema…    27     70   200 Drug      1114.  187.     90     1          180
#> 2     1 Fema…    27     70   200 Drug      1114.  187.    180     1          180
#> 3    50 Male     28     58   100 Drug       148.   10.4    90     0          180
#> 4    50 Male     28     58   100 Drug       148.   10.4   180     0          180
#> 5   100 Fema…    26     48   200 Drug       428.   35.1    90     0          180
#> 6   100 Fema…    26     48   200 Drug       428.   35.1   180     0          180
#> # ℹ 3 more variables: fit_survival <dbl>, ci_lower <dbl>, ci_upper <dbl>
```

The two engines’ confidence intervals are constructed differently under
the hood – a Wald interval on the AFT linear predictor
vs. [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)’s own
log-transform interval for the Cox model – so they’re not expected to
match exactly even when the point estimates agree closely.

[`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md)
instead returns a plain R function evaluating $`S(t \mid x)`$ at
user-specified parameters/data/times, useful for counterfactual
scenarios or plugging into simulation code:

``` r

mod_fun <- ertte_fun(mod_aft)
mod_fun(time = 90)
#>   [1] 1.926635e-01 3.662736e-01 5.444662e-01 5.444662e-01 1.157406e-01
#>   [6] 5.444662e-01 3.051399e-01 4.373140e-01 5.444662e-01 4.473582e-01
#>  [11] 1.383390e-01 4.199106e-01 3.489079e-01 4.529445e-03 4.821396e-01
#>  [16] 8.737455e-13 1.246982e-01 5.444662e-01 2.156033e-01 2.855799e-01
#>  [21] 4.771844e-01 4.028059e-01 3.121676e-01 5.444662e-01 4.333162e-01
#>  [26] 1.617855e-01 5.444662e-01 4.133145e-01 4.859960e-02 4.577562e-01
#>  [31] 5.444662e-01 4.535862e-01 3.939133e-01 6.434864e-02 5.444662e-01
#>  [36] 2.580125e-03 5.444662e-01 5.444662e-01 4.823155e-01 4.246395e-01
#>  [41] 2.851577e-01 4.079872e-01 5.444662e-01 3.695592e-01 4.716315e-01
#>  [46] 5.444662e-01 5.444662e-01 5.444662e-01 2.637058e-01 4.994996e-01
#>  [51] 5.444662e-01 5.444662e-01 3.068159e-01 7.833535e-02 4.028943e-01
#>  [56] 2.493617e-01 5.444662e-01 4.397936e-01 1.717323e-06 2.555075e-01
#>  [61] 3.534328e-01 4.404016e-04 5.444662e-01 5.444662e-01 1.804579e-11
#>  [66] 5.444662e-01 5.444662e-01 1.877829e-01 4.836900e-01 4.766591e-01
#>  [71] 3.755872e-01 4.230798e-01 5.444662e-01 5.444662e-01 5.444662e-01
#>  [76] 7.195163e-02 1.291088e-05 4.668772e-01 5.444662e-01 3.027812e-01
#>  [81] 5.444662e-01 4.522768e-01 5.444662e-01 3.254745e-01 3.675046e-01
#>  [86] 5.444662e-01 5.444662e-01 5.444662e-01 2.281600e-02 5.444662e-01
#>  [91] 4.910847e-01 5.444662e-01 8.823716e-02 2.343531e-01 5.444662e-01
#>  [96] 5.444662e-01 5.444662e-01 8.344629e-02 5.444662e-01 4.101001e-01
#> [101] 1.154691e-01 1.883223e-01 1.889679e-01 3.440030e-01 3.581465e-01
#> [106] 5.444662e-01 5.444662e-01 5.706580e-02 5.444662e-01 5.444662e-01
#> [111] 5.444662e-01 3.240047e-01 5.444662e-01 3.348738e-02 5.444662e-01
#> [116] 3.586274e-01 8.926645e-03 3.429947e-01 5.444662e-01 3.758747e-01
#> [121] 5.375257e-02 9.190491e-02 5.444662e-01 5.444662e-01 5.444662e-01
#> [126] 5.444662e-01 5.444662e-01 3.211966e-01 5.444662e-01 4.645193e-01
#> [131] 3.461114e-01 4.221832e-01 9.632716e-02 3.740232e-01 5.444662e-01
#> [136] 3.030204e-01 0.000000e+00 3.886764e-01 5.977905e-02 3.184078e-03
#> [141] 4.692625e-01 1.640226e-02 1.245934e-02 4.487511e-01 5.444662e-01
#> [146] 2.179455e-01 3.026130e-01 2.303344e-01 2.326127e-01 5.444662e-01
#> [151] 3.554163e-01 5.444662e-01 5.444662e-01 2.847818e-01 3.051431e-02
#> [156] 5.444662e-01 3.008038e-12 5.444662e-01 3.863868e-01 5.444662e-01
#> [161] 1.596356e-01 5.444662e-01 6.017923e-02 1.908546e-01 5.576975e-06
#> [166] 2.067840e-01 2.303199e-02 4.759674e-02 5.444662e-01 5.444662e-01
#> [171] 5.444662e-01 5.444662e-01 5.444662e-01 3.585787e-01 4.460558e-01
#> [176] 5.444662e-01 2.191054e-01 3.450430e-01 2.927153e-01 1.653202e-01
#> [181] 3.323767e-01 4.204893e-01 3.442207e-01 4.524353e-01 4.052131e-01
#> [186] 3.147141e-01 4.445103e-01 4.943754e-01 5.444662e-01 5.444662e-01
#> [191] 2.869081e-01 2.680967e-07 8.668939e-02 2.847002e-01 5.444662e-01
#> [196] 5.444662e-01 3.244601e-01 5.444662e-01 4.644658e-01 1.843436e-01
#> [201] 1.342335e-01 4.388579e-01 5.348463e-04 8.149026e-02 4.645789e-01
#> [206] 5.444662e-01 5.444662e-01 3.897559e-01 8.526608e-02 3.840652e-01
#> [211] 5.444662e-01 4.932070e-01 1.424974e-03 3.819995e-01 4.122573e-01
#> [216] 1.900429e-01 4.030011e-01 3.126134e-01 9.362046e-02 3.467102e-01
#> [221] 5.444662e-01 3.113831e-01 4.508048e-01 8.199120e-04 2.941038e-01
#> [226] 1.603420e-02 3.323754e-01 4.835785e-01 8.663004e-04 5.444662e-01
#> [231] 2.015472e-01 5.444662e-01 6.178393e-03 1.369293e-01 3.379028e-01
#> [236] 4.531012e-01 4.836746e-01 4.082401e-01 2.320571e-06 4.109854e-01
#> [241] 1.063449e-01 5.444662e-01 4.163054e-01 5.444662e-01 4.101713e-01
#> [246] 2.968446e-01 3.509896e-01 5.444662e-01 2.214737e-01 5.444662e-01
#> [251] 5.444662e-01 6.854909e-09 5.444662e-01 1.182199e-01 9.600126e-02
#> [256] 3.222252e-01 5.444662e-01 5.444662e-01 2.393331e-01 3.572453e-01
#> [261] 1.998976e-01 3.210723e-01 5.444662e-01 5.444662e-01 4.602965e-01
#> [266] 2.019652e-01 3.945192e-01 4.023711e-01 3.483072e-01 5.444662e-01
#> [271] 4.929390e-14 5.444662e-01 5.444662e-01 3.136765e-01 4.498156e-01
#> [276] 5.125838e-03 3.963329e-01 5.444662e-01 5.444662e-01 5.444662e-01
#> [281] 4.162068e-01 5.444662e-01 4.956018e-01 2.521316e-13 3.814192e-01
#> [286] 5.444662e-01 1.764429e-01 5.444662e-01 3.290790e-01 3.686599e-01
#> [291] 1.103549e-02 3.034758e-02 5.444662e-01 3.420684e-01 2.335831e-06
#> [296] 5.444662e-01 4.747184e-01 5.444662e-01 5.444662e-01 8.035757e-02
```

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new event
times (with realistic censoring) from either fitted model, sampling
coefficients from their asymptotic sampling distribution – the basis for
`erplots`’ visual predictive checks, if `erplots` is installed:

``` r

simulate(mod_aft, ertte_data[1:3, ], nsim = 2, seed = 8204) |>
  dplyr::select(id, aucss, sim_id, sim_time, sim_event)
#> # A tibble: 6 × 5
#>      id aucss sim_id sim_time sim_event
#>   <int> <dbl>  <int>    <dbl>     <dbl>
#> 1     1 1114.      1     16.9         1
#> 2     2  561.      1     26.8         0
#> 3     3    0       1     31.1         1
#> 4     1 1114.      2     21.3         1
#> 5     2  561.      2     10.9         1
#> 6     3    0       2    131.          1
```

## Where to go next

- [Restricted mean survival
  time](https://ertte.djnavarro.net/articles/rmst.md) and [landmark
  event probabilities](https://ertte.djnavarro.net/articles/landmark.md)
  both reduce a TTE endpoint to a single exposure-response scalar, for
  use with `erplots`’ scalar E-R plotting grammar.
- [Stepwise covariate
  modelling](https://ertte.djnavarro.net/articles/scm.md) covers
  [`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)/
  [`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
  for building up a covariate model on top of either engine.
- [`?ertte_power`](https://ertte.djnavarro.net/reference/ertte_power.md)
  covers a power-function parameterisation for continuous covariates
  (e.g. age, weight), for either engine.

## Further reading

- Kalbfleisch, J.D. and Prentice, R.L. (2002). *The Statistical Analysis
  of Failure Time Data*, 2nd edition. Wiley. A standard reference
  covering both the AFT and Cox PH model families in depth.
- Therneau, T.M. and Grambsch, P.M. (2000). *Modeling Survival Data:
  Extending the Cox Model*. Springer. Covers
  [`cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html) and the
  proportional-hazards assumption in detail.
- [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)/[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)’s
  own documentation, which
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md)/[`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  wrap directly.
