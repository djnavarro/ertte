

#' Fit an exposure-response time-to-event AFT model based on `survreg()`
#'
#' @param formula Model formula, with a `survival::Surv()` object as the
#' response, e.g. `Surv(time, event) ~ exposure`.
#' @param data Data set
#' @param dist The AFT distribution to fit, as for `survival::survreg()`.
#' Defaults to `"weibull"`. Tested and officially supported for
#' `"exponential"`, `"weibull"`, `"lognormal"`, and `"loglogistic"` --
#' see [ertte_aft_select_distribution()] for choosing among them by AIC.
#' @param ... Other arguments passed to `survival::survreg()`.
#' @returns A survreg object with extra `ertte_aft`/`ertte_model` classes
#'
#' @details The returned object has class `c("ertte_aft", "ertte_model",
#' "survreg")`: it *is* a `survreg` object, with a little extra metadata
#' attached. This means all of the usual `survreg` methods work
#' unchanged, without needing an ertte-specific equivalent -- e.g.
#' `summary()`, `coef()`, `vcov()`, `confint()`, `predict()`, `AIC()`,
#' `BIC()`, `logLik()`, and `anova()` for comparing nested models.
#' `ertte_predict()` is a separate, ertte-specific alternative to
#' `predict()` that returns survival probabilities with confidence
#' intervals in a tidy data frame; the two are complementary, not
#' competing.
#'
#' `ertte_aft()` is the AFT-specific sibling of `ertte_coxph()`, which
#' wraps `survival::coxph()` for a semi-parametric alternative. Both
#' share the `"ertte_model"` superclass, so functions that only need
#' generic operations (`update()`, `anova()`, the SCM family) work
#' unchanged across either engine; functions with AFT-specific logic
#' (e.g. `ertte_predict()`, `ertte_fun()`) dispatch via the
#' `"ertte_aft"`/`"ertte_coxph"` subclass.
#'
#' All four supported distributions are log-location-scale AFT models:
#' `log(T) = mu + scale * W`, where `mu` is the linear predictor
#' (intercept + covariates) and `W` follows a distribution that depends
#' only on `dist` (extreme-value for `"exponential"`/`"weibull"`,
#' standard normal for `"lognormal"`, standard logistic for
#' `"loglogistic"`) -- see [ertte_predict()] and `.ertte_dist_info()`.
#'
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' mod
#'
#' # other AFT distributions are also supported
#' mod_ln <- ertte_aft(Surv(time, event) ~ aucss, ertte_data, dist = "lognormal")
#' mod_ln
#'
ertte_aft <- function(formula, data, dist = "weibull", ...) {
  .ertte_check_dist(dist)
  mod <- survival::survreg(formula = formula, data = data, dist = dist, ...)
  # unlike `glm()`, `survreg()` doesn't retain the fitting data on the
  # returned object -- store it explicitly so `ertte_predict()`/
  # `ertte_fun()`/SCM (which default `newdata`/refit from `mod$data`,
  # mirroring erglm's `glm`-based equivalents) have something to fall
  # back on.
  mod$data <- data
  .as_ertte_aft(mod, dist)
}

#' Survival-probability predictions for exposure-response TTE models
#'
#' @param object An ertte model, as returned by [ertte_aft()] or
#' `ertte_coxph()` (not yet implemented)
#' @param newdata Data frame containing cases to be predicted. Defaults
#' to the data the model was fitted to.
#' @param time Numeric vector of times at which to compute survival
#' probabilities
#' @param conf_level Confidence level for the intervals
#' @param ... Passed to methods
#' @returns A tibble with one row per combination of `newdata` row and
#' `time`
#'
#' @details `ertte_predict()` is a generic, with methods for each
#' supported engine -- see [ertte_predict.ertte_aft()].
#'
#' @export
ertte_predict <- function(object, ...) {
  UseMethod("ertte_predict")
}

#' @details The `ertte_aft` method computes the linear predictor (and
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
#' A zero-row `newdata` returns a zero-row tibble with the expected
#' columns rather than erroring -- explicit here (rather than relying on
#' `predict.survreg()`'s incidental support for a zero-row `newdata`)
#' for symmetry with `ertte_predict.ertte_coxph()`, where the equivalent
#' `survival::survfit()` call genuinely does error on a zero-row
#' `newdata` (see issue #10).
#'
#' @rdname ertte_predict
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' ertte_predict(mod, ertte_data[1:5, ], time = c(30, 60, 90))
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

#' Prediction function for an exposure-response TTE model
#'
#' @param object An ertte model, as returned by [ertte_aft()] or
#' `ertte_coxph()` (not yet implemented)
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
