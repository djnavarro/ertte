
# With zero observed events (e.g. all-censored data), `coxph(model =
# TRUE)` still leaves `object$model` unset (a `survival::coxph()`
# quirk, not an ertte-introduced one) -- so the `model = TRUE` workaround
# `ertte_coxph()` otherwise relies on (see its Details) doesn't help,
# and `survival::survfit()`/`survival::basehaz()` on such a fit fail
# with a cryptic error ("'data' must be a data.frame, environment, or
# list"), from trying to re-evaluate `object$call$data` in a frame
# where that name doesn't resolve to the original fitting data. Checked
# up front in `ertte_predict.ertte_coxph()`/`ertte_fun.ertte_coxph()`/
# `.ertte_simulate_draws.ertte_coxph()` so users see an informative
# error instead. The model can still be *fit* on all-censored data
# (`ertte_coxph()` itself doesn't call this check) -- coefficients come
# back `NA` but that's a legitimate (if degenerate) result to inspect
# via `summary()`/`coef()`; only the baseline-hazard-based downstream
# methods are actually broken.
.ertte_check_coxph_nevent <- function(object) {
  if (identical(object$nevent, 0L) || identical(object$nevent, 0)) {
    rlang::abort(paste0(
      "This model has zero observed events (all-censored data), so it ",
      "has no baseline hazard to build survival predictions/simulations ",
      "from. `coef()`/`summary()` still work (coefficients are `NA` at ",
      "this degenerate fit), but `ertte_predict()`/`ertte_fun()`/",
      "`simulate()` are not supported for a zero-event `ertte_coxph` model."
    ))
  }
}

#' Cox proportional hazard modelling
#'
#' Fits a semi-parametric Cox proportional-hazards regression of
#' time-to-event data.
#'
#' @param formula Model formula specifying the regression model, 
#' e.g. `Surv(time, event) ~ exposure`.
#' @param data Data frame containing the variables of interest.
#' @param ... Other arguments passed to [survival::coxph()].
#' 
#' @return 
#' A coxph object with with additional `ertte_coxph` and `ertte_model` 
#' classes used to supply additional methods.
#'
#' @details
#' The `ertte_coxph()` function is a thin wrapper around [survival::coxph()],
#' used to build semi-parametric proportional hazard models for time-to-event 
#' data. One slight difference is that the call to `coxph()` always sets 
#' `model = TRUE`, because downstream methods supplied by the ertte package
#' require the model frame to be accessible from within the return object. 
#' Along similar lines, `ertte_coxph()` caches the input data frame within
#' the returned model object, ensuring that it remains accessible to downstream
#' tools that have access to the model object but not necessarily the original
#' data frame. 
#' 
#' Because the return value is a `coxph` object, all the usual 
#' methods for survival regression models work unchanged, without needing an 
#' ertte-specific equivalent. This includes `summary()`, `coef()`, `vcov()`, 
#' `confint()`, `predict()`, `AIC()`, `BIC()`, `logLik()`, and `anova()`.
#' Additional methods supplied via the ertte-specific classes include 
#' `ertte_predict()`, `ertte_fun()`, and `simulate()`.
#'
#' @export
#' @examples
#' mod <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' mod
#'
ertte_coxph <- function(formula, data, ...) {
  .ertte_check_response_time(formula, data)
  .ertte_check_coxph_data_size(formula, data)
  mod <- survival::coxph(formula = formula, data = data, model = TRUE, ...)
  mod$data <- data
  .as_ertte_coxph(mod)
}

