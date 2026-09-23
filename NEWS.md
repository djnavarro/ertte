# ertte 0.1.0

Initial CRAN release.

* Model fitting and prediction for exposure-response models of
  time-to-event endpoints: `ertte_aft()` (a parametric accelerated
  failure time model built on `survival::survreg()`, supporting
  exponential, Weibull, log-normal, and log-logistic distributions)
  and `ertte_coxph()` (a semi-parametric proportional-hazards model
  built on `survival::coxph()`).
* `ertte_aft_select_distribution()`, AIC-based selection across the
  supported AFT base distributions.
* `ertte_predict()` and `ertte_fun()` for survival-probability
  prediction and counterfactual curve evaluation without refitting,
  for both engines.
* `ertte_landmark()` and `ertte_rmst()`, reducing a fitted survival
  curve to a binary landmark response or a restricted mean survival
  time, with confidence intervals.
* `simulate.ertte_model()`, drawing simulated event/censoring times
  from a fitted model, with an optional per-row administrative
  `censor_time`.
* Stepwise covariate modelling (`ertte_scm_forward()`,
  `ertte_scm_backward()`, `ertte_scm_history()`), built on the
  single-term `ertte_add_term()`/`ertte_remove_term()` helpers, with
  a choice of `"p-value"`, `"aic"`, or `"bic"` selection criterion.
* `ertte_power()`, a power-function covariate transform usable as a
  model term in either engine's formula or an SCM candidate set.
* Interoperability with the companion
  [erplots](https://github.com/djnavarro/erplots) package via
  `er_predict()`/`er_simulate()`/`er_summary()`/
  `er_predict_survival()` methods, registered lazily so ertte has no
  hard dependency on erplots or on any plotting package.
* An example dataset, `ertte_data`, simulated from a Weibull AFT
  ground truth with an exposure and sex effect, administrative
  censoring, and independent dropout censoring.
