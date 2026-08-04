
#' Fit an exposure-response time-to-event Cox PH model based on `coxph()`
#'
#' @param formula Model formula, with a `survival::Surv()` object as the
#' response, e.g. `Surv(time, event) ~ exposure`.
#' @param data Data set
#' @param ... Other arguments passed to `survival::coxph()`.
#' @returns A coxph object with extra `ertte_coxph`/`ertte_model` classes
#'
#' @details The returned object has class `c("ertte_coxph", "ertte_model",
#' "coxph")`: it *is* a `coxph` object, with a little extra metadata
#' attached. This means all of the usual `coxph` methods work unchanged,
#' without needing an ertte-specific equivalent -- e.g. `summary()`,
#' `coef()`, `vcov()`, `confint()`, `predict()`, `AIC()`, `BIC()`,
#' `logLik()`, and `anova()` for comparing nested models.
#'
#' `ertte_coxph()` is the semi-parametric sibling of [ertte_aft()].
#' Unlike `ertte_aft()`, there's no `dist` argument: Cox PH doesn't
#' assume a parametric baseline hazard, so there's no distribution to
#' select. Both share the `"ertte_model"` superclass, so functions that
#' only need generic operations (`update()`, `anova()`, the SCM family)
#' work across either engine; functions with engine-specific logic (e.g.
#' `ertte_predict()`, `ertte_fun()`) dispatch via the
#' `"ertte_aft"`/`"ertte_coxph"` subclass.
#'
#' @section Not yet implemented:
#' `ertte_predict()`, `ertte_fun()`, and `simulate()` don't yet have
#' `ertte_coxph` methods -- unlike the closed-form survival function
#' available for AFT models, Cox PH prediction/simulation needs a
#' baseline hazard estimate (e.g. via `survival::survfit()`), which is
#' separate follow-up work (see AGENTS.md). Calling any of these on an
#' `ertte_coxph` object currently errors with "no applicable method".
#'
#' @export
#' @examples
#' mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
#' mod
#'
ertte_coxph <- function(formula, data, ...) {
  mod <- survival::coxph(formula = formula, data = data, ...)
  # as with `ertte_aft()`/`survreg()`, `coxph()` doesn't retain the
  # fitting data on the returned object -- store it explicitly so
  # downstream ertte functions (which default `newdata`/refit from
  # `mod$data`) have something to fall back on.
  mod$data <- data
  .as_ertte_coxph(mod)
}
