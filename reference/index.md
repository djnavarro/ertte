# Package index

## Build

Fit exposure-response time-to-event models based on survreg()/coxph()

- [`ertte_aft()`](https://ertte.djnavarro.net/reference/ertte_aft.md) :

  Fit an exposure-response time-to-event AFT model based on `survreg()`

- [`ertte_coxph()`](https://ertte.djnavarro.net/reference/ertte_coxph.md)
  :

  Fit an exposure-response time-to-event Cox PH model based on `coxph()`

- [`ertte_predict()`](https://ertte.djnavarro.net/reference/ertte_predict.md)
  : Survival-probability predictions for exposure-response TTE models

- [`ertte_landmark()`](https://ertte.djnavarro.net/reference/ertte_landmark.md)
  : Landmark event-probability predictions for exposure-response TTE
  models

- [`ertte_rmst()`](https://ertte.djnavarro.net/reference/ertte_rmst.md)
  : Restricted mean survival time predictions for exposure-response TTE
  models

- [`ertte_select_distribution()`](https://ertte.djnavarro.net/reference/ertte_select_distribution.md)
  : Select an AFT distribution by AIC

## Covariate selection

Stepwise covariate modelling for exposure-response TTE models

- [`ertte_scm_forward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
  [`ertte_scm_backward()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
  [`ertte_scm_history()`](https://ertte.djnavarro.net/reference/ertte_scm.md)
  : Stepwise covariate modelling for exposure-response TTE models
- [`ertte_add_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)
  [`ertte_remove_term()`](https://ertte.djnavarro.net/reference/ertte_term.md)
  : Add or remove a covariate term from an exposure-response TTE model
- [`ertte_power()`](https://ertte.djnavarro.net/reference/ertte_power.md)
  [`makepredictcall(`*`<ertte_power>`*`)`](https://ertte.djnavarro.net/reference/ertte_power.md)
  : Power-function covariate transform for exposure-response TTE models

## Simulate

Simulation tools for exposure-response TTE models

- [`ertte_fun()`](https://ertte.djnavarro.net/reference/ertte_fun.md) :
  Prediction function for an exposure-response TTE model
- [`simulate(`*`<ertte_model>`*`)`](https://ertte.djnavarro.net/reference/simulate.ertte_model.md)
  : Simulate from an exposure-response TTE model

## Other

Other functions and objects

- [`ertte_data`](https://ertte.djnavarro.net/reference/ertte_data.md) :
  Sample simulated data for exposure-response time-to-event models
