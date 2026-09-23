
#' Simulate from an exposure-response TTE model
#'
#' `simulate()` method for `ertte_model` objects. A single shared method
#' covers both `ertte_aft` and `ertte_coxph` fits -- there's no separate
#' `simulate.ertte_coxph()` -- with the engine-specific simulation
#' mechanics applied automatically based on the fitted object's class
#' (see Details).
#'
#' @param object An ertte model object, as returned by [ertte_aft()] or
#' [ertte_coxph()]
#' @param nsim Number of simulation replicates. Defaults to `100`.
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
#' inverse-CDF sampling, differing by engine: for
#' `ertte_aft` fits, directly from the fitted log-location-scale AFT
#' distribution (see [ertte_aft()] Details); for `ertte_coxph` fits, by
#' inverting the fitted baseline cumulative hazard ([survival::basehaz()],
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
#' that specific bias.
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



# shared generic: draws `nsim` sets of coefficients from the sampling
# distribution implied by the model's variance-covariance matrix, and
# for each draw simulates an event time per row of `newdata` via
# inverse-CDF sampling. A generic (not a single function) because the
# sampling mechanics genuinely differ by engine -- AFT samples directly
# from the fitted log-location-scale distribution
# (`.ertte_simulate_draws.ertte_aft()`), while Cox PH inverts the
# fitted baseline cumulative hazard (`.ertte_simulate_draws.ertte_coxph()`,
# in `R/ertte-coxph.R`). Used directly by `simulate.ertte_model()` (via
# `.ertte_resample()`) and `er_simulate.ertte_model()` (used by
# erplots, if installed, for TTE visual predictive checks) -- both work
# for either engine automatically via this dispatch.
#
# Administrative/observed censoring is reproduced via
# `.ertte_apply_admin_censoring()`: by default (`censor_time = NULL`),
# censored rows are capped at their *observed* exit time (the `time`
# variable in `newdata`) -- exactly correct, since that's genuinely when
# censoring happened -- while event rows are left uncensored, since
# their observed exit time is when the event happened, not their
# (unobserved) administrative censoring horizon. A genuine per-row
# administrative follow-up time can be supplied via `censor_time`
# instead, which then caps every row uniformly. `newdata` must contain
# the original response columns (`time`/`event`, named as in the
# model's `Surv()` call) -- see `.ertte_check_newdata_response()`.
.ertte_simulate_draws <- function(object, newdata, nsim = 100, seed = NULL, censor_time = NULL) {
  UseMethod(".ertte_simulate_draws")
}

.ertte_simulate_draws.ertte_aft <- function(object, newdata, nsim = 100, seed = NULL, censor_time = NULL) {
  .ertte_check_nsim(nsim)
  seed <- .ertte_pick_seed(seed)
  vars <- .ertte_check_newdata_response(object, newdata)
  censor_time <- .ertte_check_censor_time(censor_time, nrow(newdata))
  info <- .ertte_dist_info(object$ertte$type)
  scale <- object$scale
  obs_time <- newdata[[vars$time]]
  event_obs <- newdata[[vars$event]]
  withr::with_seed(
    seed = seed,
    code = {
      coef_names <- names(stats::coef(object))
      par <- mvtnorm::rmvnorm(
        n = nsim,
        mean = stats::coef(object),
        # `vcov()` also carries a row/column for `Log(scale)` (when the
        # scale is estimated jointly with the location coefficients);
        # only the location-coefficient block is needed here, since
        # `scale` itself is held fixed at its point estimate throughout
        # this package (see `ertte_aft()` Details).
        sigma = stats::vcov(object)[coef_names, coef_names, drop = FALSE]
      )
      sim <- list()
      for (ii in seq_len(nsim)) {
        dd_sim <- newdata |> dplyr::mutate(row_id = dplyr::row_number(), sim_id = ii)
        mm <- stats::model.matrix(stats::delete.response(stats::terms(object)), dd_sim)
        mu <- as.vector(mm %*% par[ii, ])
        u <- stats::runif(nrow(dd_sim))
        sim_time_raw <- exp(mu + scale * info$qbase(u))
        censored <- .ertte_apply_admin_censoring(sim_time_raw, obs_time, event_obs, censor_time)
        dd_sim$sim_time <- censored$sim_time
        dd_sim$sim_event <- censored$sim_event
        coef_draw <- stats::setNames(as.list(par[ii, ]), paste0("coef_", coef_names))
        dd_sim <- dd_sim |> dplyr::bind_cols(tibble::as_tibble(coef_draw))
        sim[[ii]] <- dd_sim
      }
    }
  )
  dplyr::bind_rows(sim)
}


# `.ertte_simulate_draws()` method for `ertte_coxph` models -- see the
# generic's documentation in `R/ertte-aft.R`. Coefficients are sampled
# from the same asymptotic normal approximation as the AFT method, but
# event times are drawn by inverting the fitted baseline cumulative
# hazard (via `.ertte_coxph_invert_basehaz()`) rather than sampling
# directly from a parametric distribution -- the baseline hazard/means
# are always taken from the fitted `object`, never recomputed for a
# sampled coefficient draw (recomputing it would need refitting the
# partial likelihood's risk sets at each draw), the same simplification
# `ertte_fun.ertte_coxph()` makes for a user-supplied `param`.
.ertte_simulate_draws.ertte_coxph <- function(object, newdata, nsim = 100, seed = NULL, censor_time = NULL) {
  .ertte_check_coxph_nevent(object)
  .ertte_check_nsim(nsim)
  seed <- .ertte_pick_seed(seed)
  vars <- .ertte_check_newdata_response(object, newdata)
  censor_time <- .ertte_check_censor_time(censor_time, nrow(newdata))
  ff <- stats::delete.response(stats::terms(object))
  means <- object$means
  bh <- survival::basehaz(object, centered = TRUE)
  obs_time <- newdata[[vars$time]]
  event_obs <- newdata[[vars$event]]
  withr::with_seed(
    seed = seed,
    code = {
      coef_names <- names(stats::coef(object))
      par <- mvtnorm::rmvnorm(
        n = nsim,
        mean = stats::coef(object),
        sigma = stats::vcov(object)[coef_names, coef_names, drop = FALSE]
      )
      sim <- list()
      for (ii in seq_len(nsim)) {
        dd_sim <- newdata |> dplyr::mutate(row_id = dplyr::row_number(), sim_id = ii)
        mm <- stats::model.matrix(ff, dd_sim)
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
        lp <- as.vector(mm %*% par[ii, ]) - as.vector(means %*% par[ii, ])
        u <- stats::runif(nrow(dd_sim))
        target_h <- -log(u) / exp(lp)
        sim_time_raw <- .ertte_coxph_invert_basehaz(bh, target_h)
        # `sim_time_raw == Inf` means the simulated draw would need to
        # survive past the fitted baseline hazard's support (the last
        # observed follow-up time across the whole cohort) to "fail" --
        # there's no information past that point either way, so treat it
        # as censored there (never as an event), matching the flat
        # extrapolation `ertte_predict.ertte_coxph()` already uses beyond
        # the observed range. Substituting `max(bh$time)` for `Inf`
        # before applying `censor_time`/`obs_time` lets a smaller cap
        # still take precedence where applicable.
        is_extrapolated <- is.infinite(sim_time_raw)
        sim_time_capped <- ifelse(is_extrapolated, max(bh$time), sim_time_raw)
        censored <- .ertte_apply_admin_censoring(sim_time_capped, obs_time, event_obs, censor_time)
        censored$sim_event[is_extrapolated] <- 0
        dd_sim$sim_time <- censored$sim_time
        dd_sim$sim_event <- censored$sim_event
        coef_draw <- stats::setNames(as.list(par[ii, ]), paste0("coef_", coef_names))
        dd_sim <- dd_sim |> dplyr::bind_cols(tibble::as_tibble(coef_draw))
        sim[[ii]] <- dd_sim
      }
    }
  )
  dplyr::bind_rows(sim)
}
