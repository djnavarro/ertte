

#' AFT regression modelling
#'
#' Fits a parametric accelerated failure time (AFT) regression model for
#' time-to-event data.
#'
#' @param formula Model formula specifying the regression model, 
#' e.g. `Surv(time, event) ~ exposure`.
#' @param data Data frame containing the variables of interest.
#' @param dist The AFT distribution type to use (default is `"weibull"`). 
#' @param ... Other arguments passed to [survival::survreg()].
#' 
#' @return 
#' A `survreg` object with additional `ertte_aft` and `ertte_model` 
#' classes used to supply additional methods.
#'
#' @details 
#' The `ertte_aft()` function is a thin wrapper around [survival::survreg()],
#' used to build parametric AFT models for time-to-event data. Unlike the original, 
#' it caches the input data frame within the returned model object, ensuring that 
#' the data set remains accessible to downstream tools that have access to the 
#' model object but not necessarily the original data frame. 
#' 
#' Like the `survreg()` function upon which it is based, `ertte_aft()` fits a
#' log-location-scale AFT model,
#'
#' \deqn{\log(t) = \mu + \sigma w}
#'
#' where \eqn{\mu} is the linear predictor, \eqn{\sigma} is the scale
#' parameter, and \eqn{w} follows a fixed base distribution determined by
#' `dist`:
#'
#' - `"weibull"`/`"exponential"`: \eqn{w} follows a standard extreme-value
#'   distribution.
#' - `"lognormal"`: \eqn{w} follows a standard normal distribution.
#' - `"loglogistic"`: \eqn{w} follows a standard logistic distribution.
#'
#' In principle `ertte_aft()` could support the full range
#' of distributions available to `survreg()`, but it is important to note that
#' the extended functionality provided by the ertte package, and hooks into other
#' exposure-response tools supplied by other packages such as erplots, have not 
#' been thoroughly tested beyond four core cases. Tested values for the `dist` 
#' argument are `weibull`, `exponential`, `lognormal`, and `loglogistic`. It is 
#' for this reason that the output object from `ertte_aft()` contains additional
#' classes: the ertte-specific tools are not necessarily valid for all `survreg`
#' objects, only those for which it has been explicitly defined.
#' 
#' Nevertheless, because the return value is genuinely a `survreg` object, albeit
#' with some additional information and metadata stored internally, all the usual 
#' methods for survival regression models work unchanged, without needing any 
#' ertte-specific equivalent. This includes `summary()`, `coef()`, `vcov()`, 
#' `confint()`, `predict()`, `AIC()`, `BIC()`, `logLik()`, and `anova()`. 
#' Additional methods supplied via the ertte-specific classes include 
#' `ertte_predict()`, `ertte_fun()`, and `simulate()`.
#'
#' @export
#' @examples
#' # fit a Weibull AFT model; because the ertte package reexports Surv(), 
#' # there is no need to load the survival package to build this model 
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' mod
#'
#' # other AFT distributions are also supported
#' mod_ln <- ertte_aft(Surv(time, event) ~ aucss, ertte_data, dist = "lognormal")
#' mod_ln
#'
ertte_aft <- function(formula, data, dist = "weibull", ...) {
  .ertte_check_dist(dist)
  mod <- survival::survreg(formula = formula, data = data, dist = dist, ...)
  # unlike `glm()`, `survreg()` doesn't retain the fitting data on the
  # returned object -- store it explicitly so `ertte_predict()`/
  # `ertte_fun()`/SCM (which default `newdata`/refit from `mod$data`,
  # mirroring erglm's `glm`-based equivalents) have something to fall
  # back on.
  mod$data <- data
  .as_ertte_aft(mod, dist)
}
