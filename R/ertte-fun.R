

#' Prediction function for an exposure-response TTE model
#'
#' Returns a closure that evaluates a fitted ertte model's survival
#' function at user-specified data, times, and (optionally) counterfactual
#' parameters, without needing to refit the model.
#'
#' @param object An ertte model, as returned by [ertte_aft()] or
#' [ertte_coxph()]
#' @param ... Passed to methods
#'
#' @returns A function with arguments `data`, `time`, and `param`, in
#' that order -- matching the argument order every other data-taking
#' entry point in the package uses (`ertte_predict()`, `ertte_landmark()`,
#' `ertte_rmst()` all take `newdata`/`data` immediately after `object`).
#' - The `data` argument should be a data frame or tibble; defaults to
#'   `object$data` (the data the model was fitted to) if not supplied.
#' - The `time` argument gives the time(s) at which to evaluate the
#'   survival function; recycled against `data`.
#' - The `param` argument should be a vector of location coefficients;
#'   defaults to `coef(object)` (the fitted coefficients) if not supplied.
#'
#' @details `ertte_fun()` is a generic, with methods for each supported
#' engine -- see [ertte_fun.ertte_aft()]. Named `ertte_fun()` for
#' consistency with the companion `erglm`/`emaxnls` packages'
#' `erglm_fun()`/`emax_fun()`, which serve the same purpose for their
#' respective model classes.
#'
#' @export
ertte_fun <- function(object, ...) {
  UseMethod("ertte_fun")
}

#' @details The `ertte_aft` method takes a fitted AFT model as input and
#' returns a function that evaluates the survival function `S(t)` at
#' user-specified parameters, data, and times (e.g. for VPCs or other
#' counterfactual simulation scenarios). The returned function checks
#' that `param` is numeric and has one entry per column of the model
#' matrix implied by `data`, erroring informatively rather than failing
#' with a cryptic "non-conformable arguments" error from matrix
#' multiplication. `scale` is always taken from the fitted `object`, not
#' from `param` (the coefficient vector from `coef()` never includes
#' it). `time` is validated the same way [ertte_predict()] validates it
#' (a numeric vector of strictly positive values) -- a non-positive
#' `time` previously returned a silent `NaN` (via `log()`) instead of
#' erroring.
#'
#' @rdname ertte_fun
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' mod_fun <- ertte_fun(mod)
#'
#' # no arguments: reproduces the fitted model's own survival predictions
#' s1 <- mod_fun(time = 60)
#'
#' # user modifies the parameters
#' par2 <- coef(mod)
#' par2["(Intercept)"] <- par2["(Intercept)"] + 1
#' s2 <- mod_fun(param = par2, time = 60)
#'
ertte_fun.ertte_aft <- function(object, ...) {
  ff <- stats::delete.response(stats::terms(object))
  info <- .ertte_dist_info(object$ertte$type)
  scale <- object$scale
  force(ff)
  function(data = NULL, time, param = NULL) {
    .ertte_check_time(time)
    if (is.null(param)) param <- stats::coef(object)
    if (is.null(data)) data <- object$data
    mm <- stats::model.matrix(ff, data)
    if (!is.numeric(param) || length(param) != ncol(mm)) {
      rlang::abort(paste0(
        "`param` must be a numeric vector of length ", ncol(mm),
        " (one entry per column of the model matrix: ",
        paste(colnames(mm), collapse = ", "), "), not length ",
        length(param), "."
      ))
    }
    mu <- as.vector(mm %*% param)
    z <- (log(time) - mu) / scale
    1 - info$pbase(z)
  }
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
#' unlike `ertte_fun.ertte_aft()`. `time` is validated the same way
#' [ertte_predict()] validates it (a numeric vector of strictly positive
#' values) -- a non-positive `time` previously returned a silent `1`
#' (as if survival were guaranteed) instead of erroring.
#'
#' @rdname ertte_fun
#' @export
#' @examples
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
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
  .ertte_check_coxph_nevent(object)
  ff <- stats::delete.response(stats::terms(object))
  means <- object$means
  bh <- survival::basehaz(object, centered = TRUE)
  force(ff)
  force(means)
  force(bh)
  function(data = NULL, time, param = NULL) {
    .ertte_check_time(time)
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

# Inverts a fitted (right-continuous, step-function) baseline
# cumulative hazard `bh` at arbitrary hazard values: the smallest
# observed `bh$time` whose cumulative hazard is at least `target_h`, or
# `Inf` if `target_h` exceeds every observed hazard value (i.e. the
# simulated event would occur after the last observed follow-up --
# left `Inf` so it gets capped/censored at the row's observed exit time
# downstream, the same administrative-censoring convention
# `.ertte_simulate_draws.ertte_aft()` uses). Used by
# `.ertte_simulate_draws.ertte_coxph()` for inverse-CDF sampling of
# event times: `S(t | x) = exp(-H0(t) * exp(lp)) = u` rearranges to
# `H0(t) = -log(u) / exp(lp)`, so inverting `H0` at that target value
# gives the simulated event time.
.ertte_coxph_invert_basehaz <- function(bh, target_h) {
  idx <- findInterval(target_h, bh$hazard) + 1L
  ifelse(idx > length(bh$hazard), Inf, bh$time[idx])
}