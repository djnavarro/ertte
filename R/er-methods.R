
# Interoperability with erplots ------------------------------------------
#
# ertte has no hard dependency on erplots (a modelling package shouldn't
# need to pull in ggplot2/patchwork). But if erplots *is* installed and
# loaded, ertte's model objects should work seamlessly with erplots'
# model-agnostic plotting API, which relies on the
# `er_predict()`/`er_simulate()`/`er_summary()` generics defined in
# erplots (see `erplots::er_model_interface`).
#
# These methods are registered lazily at load time (via `.onLoad()`
# below), so that neither erplots nor its dependencies need to be
# installed for ertte's modelling functions to work standalone.
#
# `er_predict.ertte_model()` implements phase 2's Workstream B1 --
# scalar E-R views of a TTE endpoint -- via both of the design issue's
# scalar reductions: landmark-binary (`P(event by t*)`, via
# `ertte_landmark()`, see `R/ertte-landmark.R`) and RMST (via
# `ertte_rmst()`, see `R/ertte-rmst.R`). Since erplots' `er_predict(model,
# newdata, conf_level)` contract has a fixed signature with no
# TTE-specific argument, the reduction-specific argument (`landmark_time`
# or `tau`) has to be threaded through `...`; which one is supplied
# selects which reduction runs (supplying both, or neither, errors
# informatively rather than silently picking one or doing nothing).
# `ertte_rmst()`'s own `tau` argument accepts a vector (evaluating
# multiple horizons at once), but only a single value is accepted here,
# since erplots' scalar E-R grammar expects exactly one row per `newdata`
# row. Both branches rename their reduction-specific fitted-value column
# (`fit_resp`/`fit_rmst`) to the shared `fit_resp` name erplots' plotting
# grammar expects, regardless of which reduction produced it.
#
# `er_simulate.ertte_model()` (below) mirrors this: supplying
# `landmark_time`/`tau` through `...` now also transforms its raw
# per-row `sim_time`/`sim_event` simulated draws into the `fit_resp`/
# `sim_resp` columns erplots' `er_vpc_add_simulated()` needs for a
# scalar (landmark-binary or RMST) VPC -- see `.ertte_simulate_scalar_resp()`
# for the transformation itself, and AGENTS.md's "Planned work" for the
# design decisions behind it (in particular, how genuinely ambiguous
# simulated censoring is handled).

er_predict.ertte_model <- function(model, newdata, conf_level = 0.95, ...) {
  dots <- list(...)
  has_landmark <- !is.null(dots$landmark_time)
  has_tau <- !is.null(dots$tau)

  if (has_landmark && has_tau) {
    rlang::abort(paste0(
      "`er_predict.ertte_model()` accepts only one of `landmark_time` or ",
      "`tau` (via `...`) per call -- they select different scalar E-R ",
      "reductions (landmark-binary vs restricted mean survival time)."
    ))
  }

  if (has_tau) {
    .ertte_check_single_tau(dots$tau)
    rmst <- ertte_rmst(object = model, newdata = newdata, tau = dots$tau, conf_level = conf_level)
    return(dplyr::rename(rmst, fit_resp = fit_rmst))
  }

  landmark_time <- dots$landmark_time %||% rlang::abort(paste0(
    "er_predict.ertte_model() requires either a `landmark_time` argument ",
    "(via `...`), giving the fixed time t* at which to compute P(event by ",
    "t*) -- see `ertte_landmark()` -- or a `tau` argument, giving the ",
    "restricted mean survival time horizon -- see `ertte_rmst()`."
  ))
  ertte_landmark(object = model, newdata = newdata, landmark_time = landmark_time, conf_level = conf_level)
}

er_simulate.ertte_model <- function(model, newdata, nsim = 100, seed = NULL, ...) {
  # `censor_time` isn't part of erplots' fixed `er_simulate(model,
  # newdata, nsim, seed)` contract, but can still be threaded through
  # `...` for callers that want the accurate (rather than default
  # event-rows-uncensored) simulation behaviour -- see
  # `simulate.ertte_model()`'s `censor_time` argument/Details.
  dots <- list(...)
  has_landmark <- !is.null(dots$landmark_time)
  has_tau <- !is.null(dots$tau)

  if (has_landmark && has_tau) {
    rlang::abort(paste0(
      "`er_simulate.ertte_model()` accepts only one of `landmark_time` or ",
      "`tau` (via `...`) per call -- they select different scalar E-R ",
      "reductions (landmark-binary vs restricted mean survival time)."
    ))
  }
  if (has_landmark) .ertte_check_landmark_time(dots$landmark_time)
  if (has_tau) .ertte_check_single_tau(dots$tau)

  draws <- .ertte_simulate_draws(
    object = model, newdata = newdata, nsim = nsim, seed = seed,
    censor_time = dots$censor_time
  )

  if (has_landmark || has_tau) {
    draws <- .ertte_simulate_scalar_resp(
      object = model, draws = draws,
      landmark_time = dots$landmark_time, tau = dots$tau
    )
  }
  draws
}

# Transforms `.ertte_simulate_draws()`'s raw per-(newdata row x replicate)
# `sim_time`/`sim_event` output into the `fit_resp`/`sim_resp` columns
# `erplots::er_vpc_add_simulated()` needs for a scalar (landmark-binary or
# RMST) visual predictive check (see `?erplots::er_model_interface`):
# `sim_resp` is *required* by erplots (a response-scale draw reflecting
# both parameter uncertainty and individual-level simulation noise);
# `fit_resp` is optional (a parameter-uncertainty-only draw, used for
# spaghetti-style plots), but included here too for completeness and
# parity with erglm's equivalent `er_simulate()` method.
#
# `sim_resp` is built directly from each replicate's already-censored
# `sim_time`/`sim_event` (reflecting the same administrative-censoring
# convention `.ertte_simulate_draws()` uses elsewhere -- see
# `.ertte_apply_admin_censoring()` -- which is exactly what makes the
# simulated data comparable to a real, similarly-censored observed
# study). A replicate whose simulated outcome is genuinely ambiguous
# relative to `landmark_time`/`tau` (censored strictly before it) becomes
# `NA` -- the same complete-case convention a landmark/RMST analysis
# already has to apply to genuinely censored *observed* data (e.g. the
# manual `case_when()` construction in `test-er-methods.R`'s landmark
# test), and one `er_vpc_add_simulated()`'s `mean(..., na.rm = TRUE)`
# aggregation handles correctly by simply excluding it. For RMST, the
# per-replicate individual quantity is `min(sim_time, tau)` when the
# outcome relative to `tau` is known (an event, or survival to/past
# `tau`) -- the same construction that gives `E[min(T, tau)] =
# RMST(tau)` in the population-level formalism (see the `rmst` article).
# An IPCW or pseudo-observations alternative to this complete-case
# convention was considered and rejected -- see AGENTS.md's "Planned
# work" for the full reasoning, but in short: erplots' VPC aggregation
# is an unweighted mean with nowhere for an IPCW weight to be applied,
# and pseudo-observations would need an expensive leave-one-out
# recomputation per simulated replicate.
#
# `fit_resp` reuses `ertte_fun(object)` (already implemented,
# engine-agnostic) evaluated at each replicate's own sampled coefficient
# draw -- recovered from the `coef_*` columns `.ertte_simulate_draws()`
# already attaches per replicate, so no new coefficient sampling is
# needed here. For `landmark_time` this is a single `ertte_fun()` call
# per replicate (vectorised across every `newdata` row in that replicate
# at once). For `tau` (RMST), there's no single evaluation point --
# the whole curve from 0 to `tau` needs integrating -- so
# `.ertte_rmst_fit_resp_curve()` evaluates `ertte_fun()` on a fixed grid
# and applies the composite trapezoidal rule, still vectorised across
# every `newdata` row in the replicate at once. This is a deliberately
# coarser approximation than `ertte_rmst()`'s own point estimate (an
# exact step-function sum for the Cox engine, adaptive quadrature for
# AFT) -- acceptable since `fit_resp` here is an illustrative, optional
# spaghetti-plot quantity, not something a confidence interval is built
# from.
.ertte_simulate_scalar_resp <- function(object, draws, landmark_time = NULL, tau = NULL) {
  coef_names <- names(stats::coef(object))
  coef_cols <- paste0("coef_", coef_names)
  mod_fun <- ertte_fun(object)

  draws$fit_resp <- NA_real_
  draws$sim_resp <- NA_real_

  for (ii in unique(draws$sim_id)) {
    idx <- which(draws$sim_id == ii)
    par_ii <- stats::setNames(as.numeric(draws[idx[1], coef_cols]), coef_names)
    dd <- draws[idx, , drop = FALSE]

    if (!is.null(landmark_time)) {
      draws$fit_resp[idx] <- 1 - mod_fun(param = par_ii, data = dd, time = landmark_time)
      draws$sim_resp[idx] <- dplyr::case_when(
        dd$sim_event == 1 & dd$sim_time <= landmark_time ~ 1,
        dd$sim_time > landmark_time ~ 0,
        TRUE ~ NA_real_
      )
    } else {
      draws$fit_resp[idx] <- .ertte_rmst_fit_resp_curve(mod_fun, par_ii, dd, tau)
      draws$sim_resp[idx] <- dplyr::case_when(
        dd$sim_event == 1 ~ pmin(dd$sim_time, tau),
        dd$sim_time >= tau ~ tau,
        TRUE ~ NA_real_
      )
    }
  }

  if (!is.null(landmark_time)) draws$landmark_time <- landmark_time
  if (!is.null(tau)) draws$tau <- tau
  draws
}

# Composite trapezoidal-rule integral of `mod_fun`'s survival curve from 0
# to `tau`, vectorised across every row of `data` at once (one grid point
# evaluated for the whole replicate's `newdata` rows in a single
# `mod_fun()` call, rather than one `stats::integrate()` call per row) --
# see `.ertte_simulate_scalar_resp()` for why this coarser approximation
# is an acceptable trade-off here. `S(0) = 1` always, so the grid starts
# just past 0 to avoid `log(0)` in the AFT engine's closed form.
.ertte_rmst_fit_resp_curve <- function(mod_fun, param, data, tau, n_grid = 64) {
  t_grid <- seq(0, tau, length.out = n_grid + 1)
  s_grid <- vapply(t_grid, function(tt) {
    if (tt <= 0) return(rep(1, nrow(data)))
    mod_fun(param = param, data = data, time = tt)
  }, numeric(nrow(data)))
  h <- tau / n_grid
  w <- c(1, rep(2, n_grid - 1), 1) * (h / 2)
  as.vector(s_grid %*% w)
}

er_summary.ertte_model <- function(model, conf_level = 0.95, ...) {
  coefs <- summary(model)$table
  if (is.null(coefs) || nrow(coefs) < 2) return(NULL)

  # Wald intervals, matching ertte_predict()'s approach (a normal-quantile
  # z-score applied to the standard error) rather than profile likelihood.
  z_scale <- -stats::qnorm((1 - conf_level) / 2)
  is_scale_row <- rownames(coefs) == "Log(scale)"
  coefs <- coefs[!is_scale_row, , drop = FALSE]
  estimate <- unname(coefs[, "Value"])
  std_error <- unname(coefs[, "Std. Error"])
  p_col <- grep("^p$", colnames(coefs))[1]

  coefficients <- tibble::tibble(
    term = rownames(coefs),
    estimate = estimate,
    std_error = std_error,
    statistic = if (!is.na(p_col)) unname(coefs[, p_col - 1]) else NA_real_,
    p_value = if (!is.na(p_col)) unname(coefs[, p_col]) else NA_real_,
    conf_low = estimate - z_scale * std_error,
    conf_high = estimate + z_scale * std_error,
  )

  glance <- tibble::tibble(
    n = NROW(model$y),
    df_residual = model$df.residual %||% NA_integer_,
    logLik = as.numeric(stats::logLik(model)),
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    dist = model$ertte$type,
    converged = isTRUE(model$converged) || is.null(model$converged),
  )

  p_value_slope <- if (nrow(coefficients) >= 2) coefficients$p_value[2] else NA_real_

  list(
    p_value = p_value_slope,
    coefficients = coefficients,
    glance = glance
  )
}

.onLoad <- function(libname, pkgname) {
  .s3_register("erplots::er_predict", "ertte_model", er_predict.ertte_model)
  .s3_register("erplots::er_simulate", "ertte_model", er_simulate.ertte_model)
  .s3_register("erplots::er_summary", "ertte_model", er_summary.ertte_model)
}

# Registers `method` as an S3 method for `generic` (given as
# "package::generic") and `class`, without requiring `package` to be
# installed or loaded. If `package` isn't loaded yet, registration is
# deferred until it is (via a load hook). This is the standard pattern
# used across the tidyverse for optional cross-package S3 methods (e.g.
# as implemented in `vctrs::s3_register()`); vendored here to avoid
# adding a dependency for a single small helper (same vendored copy as
# the companion `erglm` package).
.s3_register <- function(generic, class, method) {
  pieces <- strsplit(generic, "::")[[1]]
  package <- pieces[[1]]
  generic <- pieces[[2]]

  register <- function(...) {
    envir <- asNamespace(package)
    registerS3method(generic, class, method, envir = envir)
  }

  if (isNamespaceLoaded(package)) {
    register()
  }
  setHook(packageEvent(package, "onLoad"), function(...) register())

  invisible()
}
