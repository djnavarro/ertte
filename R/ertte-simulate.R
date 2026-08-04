
#' Simulate from an exposure-response TTE model
#'
#' `simulate()` method for `ertte_model` objects. Works for both
#' `ertte_aft` and `ertte_coxph` fits via a single shared method --
#' there's no separate `simulate.ertte_coxph()` -- since the
#' engine-specific simulation mechanics are dispatched internally by
#' `.ertte_simulate_draws()` (see Details).
#'
#' @param object An ertte model object, as returned by [ertte_aft()] or
#' [ertte_coxph()]
#' @param nsim Number of simulation replicates
#' @param seed Optional seed. If `NULL` (the default), one is chosen
#' automatically and reported via a message (since it determines the
#' actual simulated values returned).
#' @param newdata Data frame to simulate from. Defaults to the data the
#' model was fitted to. Must contain the original response columns
#' (`time`/`event`, as named in the model's `Surv()` call) -- see Details.
#' @param ... Unused, present for compatibility with the `simulate()`
#' generic
#'
#' @returns A tibble with one row per observation per replicate:
#' `dat_id`/`sim_id`, sampled `coef_*` columns, `sim_time` (the
#' simulated event/censoring time), and `sim_event` (1 = event, 0 =
#' censored).
#'
#' @details Coefficients are sampled from the asymptotic sampling
#' distribution implied by `vcov(object)`. Event times are then drawn by
#' inverse-CDF sampling, via the internal `.ertte_simulate_draws()` S3
#' generic, whose per-engine methods differ in exactly how: for
#' `ertte_aft` fits, directly from the fitted log-location-scale AFT
#' distribution (see [ertte_aft()] Details); for `ertte_coxph` fits, by
#' inverting the fitted baseline cumulative hazard (`survival::basehaz()`,
#' held fixed regardless of the sampled coefficient draw -- the same
#' simplification [ertte_fun.ertte_coxph()] makes for a user-supplied
#' `param`). Simulated event times are capped at each row's *observed*
#' exit time (`sim_time <- pmin(sim_time_raw, observed_time)`, with
#' `sim_event` set accordingly) to reproduce the study's observed
#' censoring/follow-up pattern -- a documented simplification, since the
#' true administrative censoring time for subjects who had an event
#' isn't otherwise available (see `.ertte_simulate_draws()`).
#'
#' @exportS3Method stats::simulate
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' sim <- simulate(mod, nsim = 20, seed = 1234)
#' sim
#'
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' sim_cox <- simulate(mod_cox, nsim = 20, seed = 1234)
#' sim_cox
#'
simulate.ertte_model <- function(object, nsim = 100, seed = NULL, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data
  .ertte_resample(object = object, newdata = newdata, nsim = nsim, seed = seed)
}

.ertte_resample <- function(object, newdata, nsim = 100, seed = NULL) {
  draws <- .ertte_simulate_draws(object = object, newdata = newdata, nsim = nsim, seed = seed)
  draws |>
    dplyr::rename(dat_id = row_id) |>
    dplyr::select(dat_id, sim_id, sim_time, sim_event, dplyr::everything())
}
