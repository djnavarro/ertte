
#' Landmark event-probability predictions for exposure-response TTE models
#'
#' Reduces the survival curve for an exposure-response time-to-even model
#' to a binary landmark event probability at a single fixed time.
#'
#' @param object An ertte model.
#' @param newdata Data frame containing cases to be predicted. Defaults
#' to the data the model was fitted to.
#' @param landmark_time A single positive number.
#' @param conf_level Confidence level for the intervals. Defaults to `.95`.
#' @param ... Passed to [ertte_predict()].
#' 
#' @return 
#' A tibble with one row per row of `newdata`. In addition to 
#' storing the original columns from `newdata`, it contains a column 
#' storing the `landmark_time` itself, a `fit_resp` column with 
#' the estimated probability of observing the event at or before the 
#' land mark time, and `ci_lower`, and `ci_upper` columns specifying 
#' the confidence interval.
#'
#' @details 
#' In some time-to-event analyses it is convenient to reduce a 
#' time-to-event endpoint to a binary landmark response: for every
#' row in the data, the landmark response describes the modelled 
#' probability that the event occurs at or before a specified 
#' `landmark_time` defined by the analyst. The landmark event
#' probability is simply one minus the survival probability at
#' that time, and is computed using `ertte_predict()`. However, 
#' unlike the typical usage of `ertte_predict()` in which a vector
#' of times is passed, `landmark_time` must be a single
#' fixed value -- a landmark is by definition evaluated at one time.
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
