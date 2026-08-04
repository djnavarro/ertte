
# With zero observed events (e.g. all-censored data), `coxph(model =
# TRUE)` still leaves `object$model` unset (a `survival::coxph()`
# quirk, not an ertte-introduced one) -- so the `model = TRUE` workaround
# `ertte_coxph()` otherwise relies on (see its Details) doesn't help,
# and `survival::survfit()`/`survival::basehaz()` on such a fit fail
# with a cryptic error ("'data' must be a data.frame, environment, or
# list"), from trying to re-evaluate `object$call$data` in a frame
# where that name doesn't resolve to the original fitting data. Checked
# up front in `ertte_predict.ertte_coxph()`/`ertte_fun.ertte_coxph()`/
# `.ertte_simulate_draws.ertte_coxph()` so users see an informative
# error instead. The model can still be *fit* on all-censored data
# (`ertte_coxph()` itself doesn't call this check) -- coefficients come
# back `NA` but that's a legitimate (if degenerate) result to inspect
# via `summary()`/`coef()`; only the baseline-hazard-based downstream
# methods are actually broken.
.ertte_check_coxph_nevent <- function(object) {
  if (identical(object$nevent, 0L) || identical(object$nevent, 0)) {
    rlang::abort(paste0(
      "This model has zero observed events (all-censored data), so it ",
      "has no baseline hazard to build survival predictions/simulations ",
      "from. `coef()`/`summary()` still work (coefficients are `NA` at ",
      "this degenerate fit), but `ertte_predict()`/`ertte_fun()`/",
      "`simulate()` are not supported for a zero-event `ertte_coxph` model."
    ))
  }
}

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
#' on the fitted baseline hazard. `simulate()` works too, via the
#' shared `simulate.ertte_model()` method -- no separate
#' `simulate.ertte_coxph()` is needed, since it dispatches internally
#' (via `.ertte_simulate_draws()`) on engine.
#'
#' @export
#' @examples
#' mod <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
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
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' ertte_predict(mod_cox, ertte_data[1:5, ], time = c(30, 60, 90))
#'
ertte_predict.ertte_coxph <- function(object, newdata = NULL, time, conf_level = .95, ...) {
  .ertte_check_coxph_nevent(object)
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

# `.ertte_simulate_draws()` method for `ertte_coxph` models -- see the
# generic's documentation in `R/ertte-core.R`. Coefficients are sampled
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
