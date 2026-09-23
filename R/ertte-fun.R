

#' Prediction function for an exposure-response TTE model
#'
#' Returns a function that evaluates a fitted ertte model's survival
#' function at user-specified data, times, and (optionally) counterfactual
#' parameters, without needing to refit the model.
#'
#' @param object An ertte model.
#' @param ... Passed to methods.
#'
#' @returns A function with arguments `data`, `time`, and `param`:
#' 
#' - The `data` argument should be a data frame or tibble; defaults to
#'   `object$data` (the data the model was fitted to) if not supplied.
#' - The `time` argument gives the time(s) at which to evaluate the
#'   survival function; recycled against `data`.
#' - The `param` argument should be a vector of location coefficients;
#'   defaults to `coef(object)` (the fitted coefficients) if not supplied.
#' 
#' The function returns a vector of survival probabilities.
#'
#' @details 
#' `ertte_fun()` is a generic function, with methods for each the supported
#' exposure-response time-to-even model classes. 
#' 
#' @section AFT models: 
#' The `ertte_aft` method takes a fitted AFT model as input and
#' returns a function that evaluates the survival function `S(t)` at
#' user-specified parameters, data, and times. Note that the `scale`
#' parameter is always taken from the fitted object, and not passed
#' via `param`. 
#' 
#' @section Cox PH models: 
#' The `ertte_coxph` method takes a fitted Cox model as input and 
#' similarly returns a function that evaluates the survival probabilities.
#' More precisely, it returns
#' 
#' \deqn{S(t|x) = S_0(t) \times \exp((x - \bar{x})' \beta)}
#' 
#' where \eqn{x} is the vector of covariates, \eqn{\bar{x}} is the mean of the covariates
#' in the fitted model, \eqn{\beta} is the vector of coefficients, and the exponentiation
#' is a matrix multiplication (i.e., \eqn{(x - \bar{x})' \beta} is the linear predictor).
#' 
#' The baseline survival function \eqn{S_0(t)} refers to the fitted  
#' baseline survival curve, calulated by calling `survival::basehaz()`
#' with `centered = TRUE` and held constant beyond the last observed 
#' time, consistent with the approach adopted by `ertte_predict()` for
#' Cox proportional hazard models. 
#' 
#' Analogous to the AFT model case, `param` is only used to vary the
#' linear predictor: the baseline hazard is always taken from the 
#' fitted model, and not recomputed for a hypothetical `param` value 
#' as that would require refitting the model. 
#' 
#' @rdname ertte_fun
#' @export
#' 
#' @examples
#' # fit an AFT model
#' mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' mod_aft_fun <- ertte_fun(mod_aft)
#'
#' # when called with no arguments it reproduces the fitted 
#' # model's own survival predictions
#' s1 <- mod_aft_fun(time = 60)
#' s1[1:5]
#'
#' # when called with an explicit parameter vector, the model
#' # is evaluated at those modified parameters
#' par_new <- coef(mod_aft)
#' par_new["(Intercept)"] <- par_new["(Intercept)"] + 1
#' s2 <- mod_aft_fun(param = par_new, time = 60)
#' s2[1:5]
#'
#' # the same logic applies for Cox regression models 
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' mod_cox_fun <- ertte_fun(mod_cox)
#' s3 <- mod_cox_fun(time = 60)
#' s3[1:5]
#' 
ertte_fun <- function(object, ...) {
  UseMethod("ertte_fun")
}

#' @name ertte_fun
#' @export
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

#' @rdname ertte_fun
#' @export
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