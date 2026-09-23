# Changelog

## ertte 0.1.0

Initial CRAN release.

- Model fitting and prediction for exposure-response models of
  time-to-event endpoints:
  [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) (a
  parametric accelerated failure time model built on
  [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html),
  supporting exponential, Weibull, log-normal, and log-logistic
  distributions) and
  [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  (a semi-parametric proportional-hazards model built on
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)).
- [`ertte_aft_select_distribution()`](https://ertte.djnavarro.net/reference/ertte_aft_select_distribution.md),
  AIC-based selection across the supported AFT base distributions.
- [`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
  and
  [`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md)
  for survival-probability prediction and counterfactual curve
  evaluation without refitting, for both engines.
- [`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)
  and
  [`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md),
  reducing a fitted survival curve to a binary landmark response or a
  restricted mean survival time, with confidence intervals.
- [`simulate.ertte_model()`](https://ertte.djnavarro.net/reference/simulate.ertte_model.md),
  drawing simulated event/censoring times from a fitted model, with an
  optional per-row administrative `censor_time`.
- Stepwise covariate modelling
  ([`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md),
  [`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md),
  [`ertte_scm_history()`](https://ertte.djnavarro.net/reference/ertte_scm.md)),
  built on the single-term
  [`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)/[`ertte_remove_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)
  helpers, with a choice of `"p-value"`, `"aic"`, or `"bic"` selection
  criterion.
- [`ertte_power()`](https://ertte.djnavarro.net/reference/ertte_power.md),
  a power-function covariate transform usable as a model term in either
  engine’s formula or an SCM candidate set.
- Interoperability with the companion
  [erplots](https://github.com/djnavarro/erplots) package via
  `er_predict()`/`er_simulate()`/`er_summary()`/ `er_predict_survival()`
  methods, registered lazily so ertte has no hard dependency on erplots
  or on any plotting package.
- An example dataset, `ertte_data`, simulated from a Weibull AFT ground
  truth with an exposure and sex effect, administrative censoring, and
  independent dropout censoring.
