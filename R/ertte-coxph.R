
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
#' `ertte_predict()` and `ertte_fun()` have `ertte_coxph` methods (see
#' [ertte_predict.ertte_coxph()]/[ertte_fun.ertte_coxph()]), both built
#' on the fitted baseline hazard. `simulate()` doesn't yet have an
#' `ertte_coxph` method -- see AGENTS.md. Calling it on an `ertte_coxph`
#' object currently errors with "no applicable method".
#'
#' @export
#' @examples
#' mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
#' mod
#'
ertte_coxph <- function(formula, data, ...) {
  # `model = TRUE` (retaining the model frame on the fitted object) is
  # required, not just a nice-to-have: `survival::survfit()` (which
  # `ertte_predict.ertte_coxph()` relies on) otherwise tries to
  # reconstruct the model frame by re-evaluating `object$call$data` --
  # which fails here for the same reason `stats::update()` fails on
  # ertte_aft()`/`ertte_coxph()` fits (see `.ertte_refit()` in
  # `R/ertte-scm.R`): the captured call refers to this function's own
  # local `formula`/`data` bindings, not anything visible in
  # `survfit()`'s caller's frame. Storing the model frame directly
  # sidesteps that.
  mod <- survival::coxph(formula = formula, data = data, model = TRUE, ...)
  # as with `ertte_aft()`/`survreg()`, `coxph()` doesn't retain the
  # fitting data on the returned object -- store it explicitly so
  # downstream ertte functions (which default `newdata`/refit from
  # `mod$data`) have something to fall back on.
  mod$data <- data
  .as_ertte_coxph(mod)
}

#' @details The `ertte_coxph` method delegates to
#' `survival::survfit(object, newdata, conf.int = conf_level)`, which
#' computes a per-row survival curve `S(t | x) = S0(t)^exp(lp(x) -
#' lp(xbar))` from the fitted baseline hazard (Breslow or Efron,
#' matching `object$method`) and the linear predictor, then evaluates it
#' at `time` via `summary(..., extend = TRUE)` -- `extend = TRUE` allows
#' `time` to exceed the last observed follow-up time, holding survival
#' constant beyond it (the usual step-function extrapolation) rather
#' than erroring. Confidence intervals come from `survfit()`'s own
#' `conf.type = "log"` transform (Wald on `log(-log(S))`), which is
#' better suited to a probability bounded in `[0, 1]` than the plain
#' Wald interval `ertte_predict.ertte_aft()` uses on the linear
#' predictor -- the two methods' intervals are not directly comparable
#' as a result, which is expected given the different model structures.
#'
#' @rdname ertte_predict
#' @export
#' @examples
#' mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
#' ertte_predict(mod_cox, ertte_data[1:5, ], time = c(30, 60, 90))
#'
ertte_predict.ertte_coxph <- function(object, newdata = NULL, time, conf_level = .95, ...) {
  .ertte_check_conf_level(conf_level)
  if (is.null(newdata)) newdata <- object$data
  if (!is.numeric(time) || length(time) == 0L || anyNA(time) || any(time <= 0)) {
    rlang::abort("`time` must be a numeric vector of strictly positive values.")
  }
  n <- nrow(newdata)
  k <- length(time)

  sf <- survival::survfit(object, newdata = newdata, conf.int = conf_level, se.fit = TRUE)
  ss <- summary(sf, times = time, extend = TRUE)
  # `summary()`'s `$surv`/`$lower`/`$upper` are `[k x n]` matrices when
  # `newdata` has more than one row, but drop to a plain length-`k`
  # vector when it has exactly one -- normalise both to a `[k x n]`
  # matrix so the flattening below doesn't need a special case.
  as_km <- function(x) matrix(x, nrow = k, ncol = n)

  rep_rows <- rep(seq_len(n), each = k)
  time_rep <- rep(time, times = n)

  out <- newdata[rep_rows, , drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      time = unname(time_rep),
      fit_survival = unname(as.vector(as_km(ss$surv))),
      ci_lower = unname(as.vector(as_km(ss$lower))),
      ci_upper = unname(as.vector(as_km(ss$upper))),
    )
  return(out)
}

# Evaluates a fitted (right-continuous, step-function) baseline
# cumulative hazard `bh` (as returned by `survival::basehaz()`, a data
# frame with `time`/`hazard` columns sorted ascending by `time`) at
# arbitrary times, held constant beyond the last observed time --
# matching `ertte_predict.ertte_coxph()`'s `extend = TRUE` behaviour so
# `ertte_fun.ertte_coxph()`'s counterfactual evaluation is consistent
# with the model's own predictions.
.ertte_coxph_basehaz_at <- function(bh, time) {
  idx <- findInterval(time, bh$time)
  ifelse(idx == 0, 0, bh$hazard[idx])
}

#' @details The `ertte_coxph` method returns a function that evaluates
#' `S(t | x) = S0(t)^exp((x - xbar)'param)`, where `S0(t)` is the fitted
#' baseline survival curve (via `survival::basehaz(object, centered =
#' TRUE)`, held constant beyond the last observed time, matching
#' `ertte_predict.ertte_coxph()`) and `xbar` is `object$means` (the
#' covariate means `coxph()` centers the partial likelihood on when
#' fitting -- centering matters here because `basehaz()`'s baseline is
#' defined relative to it, not to `x = 0`). As with
#' `ertte_fun.ertte_aft()`, `param` only varies the linear predictor:
#' the baseline hazard is always taken from the fitted `object`, never
#' recomputed for a hypothetical `param` (that would need refitting the
#' partial likelihood's risk sets) -- matching the level of
#' approximation used elsewhere in this package (e.g. `scale` for AFT
#' models is likewise held fixed). Since Cox models have no intercept,
#' `param` has one entry per covariate with no `"(Intercept)"` column,
#' unlike `ertte_fun.ertte_aft()`.
#'
#' @rdname ertte_fun
#' @export
#' @examples
#' mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
#' mod_cox_fun <- ertte_fun(mod_cox)
#'
#' # no arguments: reproduces the fitted model's own survival predictions
#' s1 <- mod_cox_fun(time = 60)
#'
#' # user modifies the parameters
#' par2 <- coef(mod_cox)
#' par2["aucss"] <- par2["aucss"] * 1.5
#' s2 <- mod_cox_fun(param = par2, time = 60)
#'
ertte_fun.ertte_coxph <- function(object, ...) {
  ff <- stats::delete.response(stats::terms(object))
  means <- object$means
  bh <- survival::basehaz(object, centered = TRUE)
  force(ff)
  force(means)
  force(bh)
  function(param = NULL, data = NULL, time) {
    if (is.null(param)) param <- stats::coef(object)
    if (is.null(data)) data <- object$data
    mm <- stats::model.matrix(ff, data)
    # `coxph()` models have no intercept (it cancels out of the partial
    # likelihood and is absorbed into the baseline hazard), but
    # `model.matrix()` on `ff` adds one anyway since the underlying
    # `terms()` object doesn't record that -- drop it so `ncol(mm)`
    # matches `length(coef(object))`.
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    if (!is.numeric(param) || length(param) != ncol(mm)) {
      rlang::abort(paste0(
        "`param` must be a numeric vector of length ", ncol(mm),
        " (one entry per column of the model matrix: ",
        paste(colnames(mm), collapse = ", "), "), not length ",
        length(param), "."
      ))
    }
    lp <- as.vector(mm %*% param) - as.vector(means %*% param)
    haz <- .ertte_coxph_basehaz_at(bh, time)
    exp(-haz * exp(lp))
  }
}
