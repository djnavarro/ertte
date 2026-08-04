
#' Landmark event-probability predictions for exposure-response TTE models
#'
#' @param object An ertte model, as returned by [ertte_aft()] or
#' [ertte_coxph()].
#' @param newdata Data frame containing cases to be predicted. Defaults
#' to the data the model was fitted to.
#' @param landmark_time A single, strictly positive number: the fixed
#' time `t*` at which to compute `P(event by t*)`.
#' @param conf_level Confidence level for the intervals.
#' @param ... Passed to [ertte_predict()].
#' @returns A tibble with one row per row of `newdata`, plus
#' `landmark_time`, `fit_resp` (the estimated `P(event by t*)`),
#' `ci_lower`, and `ci_upper`.
#'
#' @details Reduces a TTE endpoint to a binary landmark response --
#' "did the event happen by a fixed time t*" -- turning it into an
#' ordinary scalar exposure-response value that erplots' existing
#' `er_plot()`/`er_vpc()` grammars can visualise with no new plotting
#' code (see the package's design issue, Workstream B1:
#' <https://github.com/djnavarro/ertte/issues/1>). `P(event by t*) = 1
#' - S(t*)`, computed by calling [ertte_predict()] at `time =
#' landmark_time` and transforming its survival-probability output.
#' Since that's a decreasing monotonic transform, the confidence
#' interval bounds swap (the upper bound on survival becomes the lower
#' bound on event probability, and vice versa) but need no
#' recomputation of their own: whatever validity `ertte_predict()`'s
#' interval has for a given engine -- a Wald interval on the AFT
#' method's linear predictor, or `survival::survfit()`'s own
#' `conf.type = "log"` interval for the Cox PH method -- carries
#' through unchanged.
#'
#' `ertte_landmark()` is a single function, not a generic -- unlike
#' [ertte_predict()]/[ertte_fun()], it needs no engine-specific logic
#' of its own: it delegates entirely to `ertte_predict()`, which
#' already dispatches on the `ertte_aft`/`ertte_coxph` subclass. This
#' also means all of `ertte_predict()`'s existing edge-case handling
#' (e.g. the all-censored-Cox guard, single-stratum `NA` propagation)
#' is inherited unchanged.
#'
#' Unlike `ertte_predict()`'s `time` argument (a vector, evaluated at
#' potentially many times per row), `landmark_time` must be a single
#' fixed value -- a landmark is by definition evaluated at one time.
#'
#' Restricted mean survival time (RMST), the other scalar E-R
#' reduction the design issue mentions, is implemented separately as
#' [ertte_rmst()] -- unlike `ertte_landmark()`, it's a genuine generic
#' rather than a thin wrapper around `ertte_predict()`, since an area
#' under the curve needs the whole survival curve, not a single time
#' point.
#'
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' ertte_landmark(mod, ertte_data[1:5, ], landmark_time = 180)
#'
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' ertte_landmark(mod_cox, ertte_data[1:5, ], landmark_time = 180)
#'
ertte_landmark <- function(object, newdata = NULL, landmark_time, conf_level = .95, ...) {
  .ertte_check_landmark_time(landmark_time)
  pred <- ertte_predict(object, newdata = newdata, time = landmark_time, conf_level = conf_level, ...)
  pred |>
    dplyr::mutate(
      landmark_time = time,
      fit_resp = 1 - fit_survival,
      # a decreasing transform swaps which survival-scale bound
      # becomes which event-probability-scale bound
      new_ci_lower = 1 - ci_upper,
      new_ci_upper = 1 - ci_lower,
    ) |>
    dplyr::select(-time, -fit_survival, -ci_lower, -ci_upper) |>
    dplyr::rename(ci_lower = new_ci_lower, ci_upper = new_ci_upper)
}
