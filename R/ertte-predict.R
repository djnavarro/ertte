
#' Survival-probability predictions for exposure-response TTE models
#'
#' Computes fitted survival probabilities `S(t)` and confidence intervals
#' from a fitted ertte model, for one or more rows of `newdata` at one or
#' more `time` values.
#'
#' @param object An ertte model
#' @param newdata Data frame containing cases to be predicted. Defaults
#' to the data the model was fitted to.
#' @param time Numeric vector of times at which to compute survival
#' probabilities
#' @param conf_level Confidence level for the intervals. Defaults to `.95`.
#' @param ... Passed to methods
#' @returns A tibble with one row per combination of `newdata` row and
#' `time`
#'
#' @details `ertte_predict()` is a generic function, with two supplied 
#' methods, one for parametric AFT exposure-response models (i.e., 
#' `ertte_aft` classed objects) and another for Cox proportional hazard
#' exposure-response models (i.e., objects with class `ertte_coxph`). 
#'
#' The `ertte_aft` method computes the linear predictor (and
#' its standard error) via `predict(object, newdata, type = "linear",
#' se.fit = TRUE)`, then converts to a survival probability `S(t) = 1 -
#' F((log(t) - mu) / scale)`, where `F` is the base distribution's CDF
#' implied by `object`'s `dist` (see [ertte_aft()] Details). Confidence
#' intervals are Wald intervals on `mu` (a `qnorm()` z-score times the
#' standard error), back-transformed the same way -- parameter
#' uncertainty in `scale` is not propagated, matching the level of
#' approximation used throughout this package (e.g. `erglm_predict()`'s
#' equivalent in the companion `erglm` package). `conf_level` must be a
#' single number between 0 and 1 (inclusive); other values error rather
#' than silently producing a reversed or `NaN` interval.
#' 
#' The `ertte_coxph` method delegates to
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
#' # predictions for an AFT model
#' mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' ertte_predict(mod_aft, ertte_data[1:5, ], time = c(30, 60, 90))
#' 
#' # predictions for a Cox model
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' ertte_predict(mod_cox, ertte_data[1:5, ], time = c(30, 60, 90))
#' 
ertte_predict <- function(object, ...) {
  UseMethod("ertte_predict")
}



# ertte_predict method for aft models ------

#' @export
#' @name ertte_predict
#' 
ertte_predict.ertte_aft <- function(object, newdata = NULL, time, conf_level = .95, ...) {
  .ertte_check_conf_level(conf_level)
  if (is.null(newdata)) newdata <- object$data
  .ertte_check_time(time)
  if (nrow(newdata) == 0L) {
    return(
      newdata |>
        tibble::as_tibble() |>
        dplyr::mutate(
          time = numeric(0),
          fit_survival = numeric(0),
          ci_lower = numeric(0),
          ci_upper = numeric(0)
        )
    )
  }
  info <- .ertte_dist_info(object$ertte$type)
  z_scale <- -stats::qnorm((1 - conf_level) / 2)
  scale <- object$scale

  lp <- stats::predict(object, newdata, type = "linear", se.fit = TRUE)
  n <- nrow(newdata)
  k <- length(time)
  rep_rows <- rep(seq_len(n), each = k)

  time_rep <- rep(time, times = n)
  mu_rep <- rep(lp$fit, each = k)
  se_mu_rep <- rep(lp$se.fit, each = k)
  z <- (log(time_rep) - mu_rep) / scale

  out <- newdata[rep_rows, , drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      time = unname(time_rep),
      fit_survival = unname(1 - info$pbase(z)),
      ci_lower = unname(1 - info$pbase((log(time_rep) - (mu_rep - z_scale * se_mu_rep)) / scale)),
      ci_upper = unname(1 - info$pbase((log(time_rep) - (mu_rep + z_scale * se_mu_rep)) / scale)),
    )
  return(out)
}



# ertte_predict method for coxph models ------

#' @export
#' @name ertte_predict
#' 
ertte_predict.ertte_coxph <- function(object, newdata = NULL, time, conf_level = .95, ...) {
  .ertte_check_coxph_nevent(object)
  .ertte_check_conf_level(conf_level)
  if (is.null(newdata)) newdata <- object$data
  .ertte_check_time(time)
  if (nrow(newdata) == 0L) {
    return(
      newdata |>
        tibble::as_tibble() |>
        dplyr::mutate(
          time = numeric(0),
          fit_survival = numeric(0),
          ci_lower = numeric(0),
          ci_upper = numeric(0)
        )
    )
  }
  n <- nrow(newdata)
  k <- length(time)

  # `survfit()`'s own `conf.int` argument rejects exactly 0/1 -- pass a
  # harmless placeholder in that case and construct the boundary
  # interval manually below, since `$surv` (the point estimate) doesn't
  # depend on `conf.int` at all (confirmed empirically: only `$lower`/
  # `$upper` do).
  sf <- survival::survfit(
    object, newdata = newdata,
    conf.int = if (conf_level %in% c(0, 1)) 0.95 else conf_level,
    se.fit = TRUE
  )
  ss <- summary(sf, times = time, extend = TRUE)
  # `summary()`'s `$surv`/`$lower`/`$upper` are `[k x n]` matrices when
  # `newdata` has more than one row, but drop to a plain length-`k`
  # vector when it has exactly one -- normalise both to a `[k x n]`
  # matrix so the flattening below doesn't need a special case.
  as_km <- function(x) matrix(x, nrow = k, ncol = n)
  fit_mat <- as_km(ss$surv)
  if (conf_level == 0) {
    lower_mat <- fit_mat
    upper_mat <- fit_mat
  } else if (conf_level == 1) {
    lower_mat <- matrix(0, nrow = k, ncol = n)
    upper_mat <- matrix(1, nrow = k, ncol = n)
  } else {
    lower_mat <- as_km(ss$lower)
    upper_mat <- as_km(ss$upper)
  }

  rep_rows <- rep(seq_len(n), each = k)
  time_rep <- rep(time, times = n)

  out <- newdata[rep_rows, , drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      time = unname(time_rep),
      fit_survival = unname(as.vector(fit_mat)),
      ci_lower = unname(as.vector(lower_mat)),
      ci_upper = unname(as.vector(upper_mat)),
    )
  return(out)
}
