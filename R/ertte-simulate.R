
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
#' @param censor_time Optional administrative/maximum-follow-up time(s)
#' to cap simulated event times at, applied uniformly to every row
#' regardless of whether that row observed an event. Either `NULL` (the
#' default -- see Details for the fallback behaviour), a single number
#' (recycled across all rows of `newdata`), or a numeric vector of
#' length `nrow(newdata)` giving each row's own administrative follow-up
#' time.
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
#' `param`).
#'
#' The resulting *raw* (uncensored) simulated event time is then censored
#' by `censor_time` if supplied (`sim_time <- pmin(sim_time_raw,
#' censor_time)`, with `sim_event` set accordingly) -- this is the
#' accurate case, whenever a genuine per-row (or study-wide constant)
#' administrative follow-up time is known, since it caps every row (event
#' or censored) against its true censoring horizon.
#'
#' Absent a supplied `censor_time` (the default, `NULL`), row's own
#' observed `event` status determines the fallback: rows that were
#' *censored* in `newdata` have their observed exit time used as the cap
#' (`sim_time <- pmin(sim_time_raw, observed_time)`), since that
#' observed exit time genuinely is when censoring happened -- an exact
#' match, not an approximation. Rows that had an observed *event*,
#' however, are left **uncensored** in the simulation: their observed
#' exit time is when the event actually happened, not their
#' administrative censoring horizon (which was necessarily later, and
#' typically isn't recorded once an event has occurred) -- capping
#' simulated draws there would leak the observed event day into the
#' simulation and bias a simulated-vs-observed comparison (e.g. a visual
#' predictive check) toward looking more similar than the fitted model
#' actually implies. This remains an approximation for event rows (no
#' censoring is applied at all, absent better information), but avoids
#' that specific bias -- see `.ertte_apply_admin_censoring()`.
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
#' # a genuine per-row administrative censoring time -- ertte_data's
#' # `admin_censor` column records the fixed 180-day study cutoff used
#' # to generate it, known regardless of whether a subject had an event
#' sim_admin <- simulate(mod, nsim = 20, seed = 1234, censor_time = ertte_data$admin_censor)
#' sim_admin
#'
simulate.ertte_model <- function(object, nsim = 100, seed = NULL, newdata = NULL, censor_time = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data
  .ertte_resample(object = object, newdata = newdata, nsim = nsim, seed = seed, censor_time = censor_time)
}

.ertte_resample <- function(object, newdata, nsim = 100, seed = NULL, censor_time = NULL) {
  draws <- .ertte_simulate_draws(
    object = object, newdata = newdata, nsim = nsim, seed = seed, censor_time = censor_time
  )
  draws |>
    dplyr::rename(dat_id = row_id) |>
    dplyr::select(dat_id, sim_id, sim_time, sim_event, dplyr::everything())
}
