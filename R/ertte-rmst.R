
# Reimplements (rather than depends on) `survival:::survmean()`'s
# rectangle/tail-weighted-sum construction for a single survival curve's
# restricted mean, generalised to use a *profile-specific* variance
# increment instead of `survmean()`'s population-level Greenwood term.
#
# `survmean()` (used internally by `print.survfit(fit, rmean = tau)`)
# estimates the variance-of-mean via `hh <- n.event / (n.risk * (n.risk -
# n.event))` -- the classic Greenwood increment for a *plain* KM curve.
# Applied to a `coxph`-derived, per-covariate-profile curve (as returned by
# `survival::survfit(coxph_object, newdata = ...)`), `n.risk`/`n.event` are
# shared across every profile (they describe the baseline cohort's risk
# sets, not a specific covariate value) -- so that term captures baseline-
# hazard sampling variability but not coefficient uncertainty, and can
# badly *understate* the true variance for a profile far from the mean
# covariate values (confirmed empirically: for an extreme covariate
# profile, `survmean()`'s SE was ~15x smaller than a nonparametric
# bootstrap SE for the same quantity).
#
# The fix keeps the same rectangle/tail-weighted-sum machinery, but swaps
# the variance-increment source for `diff(std.err^2)`: `survfit.coxph()`'s
# own `std.err` field is already, per profile, the correctly-combined
# (coefficient + baseline-hazard) standard error of the cumulative hazard
# H(t|x) (confirmed empirically: `sf$logse == TRUE` and `sf$cumhaz ==
# -log(sf$surv)` exactly, so `std.err(t)` genuinely estimates SE[H(t|x)]).
# Treating its increments as if they accumulate independently across jump
# times (the same simplifying assumption the classic Greenwood formula
# makes) is an approximation -- the coefficient-uncertainty component of
# H(t|x) is really a single random direction shared across all t, not a
# sum of independent per-jump increments -- but cross-validated against a
# 300-replicate nonparametric bootstrap across two very different
# covariate profiles, this approximation tracked the bootstrap SE
# substantially more closely than either `survmean()`'s naive population
# term or a cheaper approximation that ignores baseline-hazard uncertainty
# entirely (holding it fixed and only redrawing coefficients).
#
# `time`/`surv`/`std.err` are one curve's worth of `survfit()` output
# (i.e. already selected out of the `[time x profile]` matrix `survfit()`
# returns when `newdata` has more than one row); `start.time` is `sf$t0`
# if present, else `min(0, sf$time)` (matching `survmean()`); `end.time`
# is the RMST horizon `tau`.
.ertte_rmst_pfun_delta <- function(time, surv, std.err, start.time, end.time) {
  varH <- std.err^2
  dVarH <- diff(c(0, varH)) # increment of Var[H(t|x)] at each jump time

  keep <- which(time <= end.time)
  if (length(keep) == 0) {
    temptime <- end.time
    tempsurv <- 1
    hh <- 0
  } else {
    temptime <- c(time[keep], end.time)
    tempsurv <- c(surv[keep], surv[max(keep)])
    hh <- c(dVarH[keep], 0)
  }
  n <- length(temptime)
  delta <- diff(c(start.time, temptime))
  rectangles <- delta * c(1, tempsurv[-n])
  varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1])
  mean <- sum(rectangles)
  c(rmean = mean, se_rmean = sqrt(varmean))
}

#' Restricted mean survival time predictions for exposure-response TTE models
#'
#' @param object An ertte model, as returned by [ertte_aft()] or
#' [ertte_coxph()].
#' @param newdata Data frame containing cases to be predicted. Defaults
#' to the data the model was fitted to.
#' @param tau Numeric vector of restriction horizons at which to compute
#' `RMST(tau) = integral of S(t) from 0 to tau`.
#' @param conf_level Confidence level for the intervals.
#' @param ... Passed to methods.
#' @returns A tibble with one row per combination of `newdata` row and
#' `tau`, plus `fit_rmst`, `ci_lower`, and `ci_upper`.
#'
#' @details `ertte_rmst()` reduces a TTE endpoint to a scalar
#' exposure-response value -- restricted mean survival time, the area
#' under the survival curve up to a fixed horizon `tau` -- the other
#' scalar reduction the package's design issue mentions alongside
#' landmark-binary (see [ertte_landmark()]). Unlike `ertte_landmark()`,
#' `ertte_rmst()` is a generic (not a thin wrapper around
#' `ertte_predict()`): computing an area under the curve genuinely needs
#' engine-specific logic (the whole curve, not a single time point), so
#' there are per-engine methods -- see [ertte_rmst.ertte_aft()] and
#' [ertte_rmst.ertte_coxph()].
#'
#' Confidence intervals are symmetric Wald intervals on the RMST scale
#' (`fit_rmst +/- z * se_rmst`) for both engines. This is not
#' automatically bounded to `[0, tau]` the way `ertte_predict()`'s
#' survival-probability intervals are bounded to `[0, 1]` by construction
#' (their back-transform through a CDF keeps them there) -- an
#' unclipped Wald interval on RMST can, in principle, dip below 0 or
#' exceed `tau` for small samples or near-boundary cases. `conf_level`
#' must be a single number between 0 and 1 (inclusive); other values
#' error rather than silently producing a reversed or `NaN` interval.
#'
#' A zero-row `newdata` returns a zero-row tibble with the expected
#' columns for both engines, rather than erroring (see issue #10 for
#' why this needed an explicit guard on the `ertte_coxph` side).
#'
#' @export
ertte_rmst <- function(object, newdata = NULL, tau, conf_level = .95, ...) {
  UseMethod("ertte_rmst")
}

#' @details The `ertte_aft` method computes `fit_rmst` by numerically
#' integrating the closed-form survival function `S(t|x) = 1 - F((log(t)
#' - mu) / scale)` from 0 to `tau` via `stats::integrate()`, where `mu`
#' (and its standard error) comes from `predict(object, newdata, type =
#' "linear", se.fit = TRUE)`, matching [ertte_predict.ertte_aft()]. The
#' standard error is an analytic delta method that differentiates under
#' the integral sign: `d/dmu RMST(tau|x) = integral of dbase(z) / scale
#' from 0 to tau`, where `dbase` is the base distribution's density (see
#' `.ertte_dist_info()`) -- propagating only `Var(mu)`, not `Var(scale)`,
#' the same simplification `ertte_predict.ertte_aft()` already makes for
#' its own confidence intervals.
#'
#' @rdname ertte_rmst
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' ertte_rmst(mod, ertte_data[1:5, ], tau = c(60, 90))
#'
ertte_rmst.ertte_aft <- function(object, newdata = NULL, tau, conf_level = .95, ...) {
  .ertte_check_conf_level(conf_level)
  if (is.null(newdata)) newdata <- object$data
  if (!is.numeric(tau) || length(tau) == 0L || anyNA(tau) || any(tau <= 0)) {
    rlang::abort("`tau` must be a numeric vector of strictly positive values.")
  }
  if (nrow(newdata) == 0L) {
    return(
      newdata |>
        tibble::as_tibble() |>
        dplyr::mutate(
          tau = numeric(0),
          fit_rmst = numeric(0),
          ci_lower = numeric(0),
          ci_upper = numeric(0)
        )
    )
  }
  info <- .ertte_dist_info(object$ertte$type)
  scale <- object$scale
  z_scale <- -stats::qnorm((1 - conf_level) / 2)

  lp <- stats::predict(object, newdata, type = "linear", se.fit = TRUE)
  n <- nrow(newdata)
  k <- length(tau)
  rep_rows <- rep(seq_len(n), each = k)
  tau_rep <- rep(tau, times = n)
  mu_rep <- rep(lp$fit, each = k)
  se_mu_rep <- rep(lp$se.fit, each = k)

  s_fun <- function(t, mu) 1 - info$pbase((log(t) - mu) / scale)
  grad_fun <- function(t, mu) info$dbase((log(t) - mu) / scale) / scale

  # named distinctly from the `fit_rmst`/`ci_lower`/`ci_upper` columns the
  # `mutate()` below creates -- `newdata` may already carry columns of
  # those exact names (e.g. when a previous `ertte_rmst()` call's own
  # output is reused as `newdata`, as `erplots::er_plot()`'s model-curve
  # grid naturally does), and `dplyr::mutate()`'s data mask would resolve
  # a bare `fit_rmst` reference against that pre-existing column instead
  # of this freshly computed vector, silently passing through stale
  # values. See erplots#12.
  fit_rmst_val <- se_rmst_val <- numeric(length(tau_rep))
  for (i in seq_along(tau_rep)) {
    fit_rmst_val[i] <- stats::integrate(s_fun, 0, tau_rep[i], mu = mu_rep[i])$value
    grad <- stats::integrate(grad_fun, 0, tau_rep[i], mu = mu_rep[i])$value
    se_rmst_val[i] <- abs(grad) * se_mu_rep[i]
  }

  newdata[rep_rows, , drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      tau = unname(tau_rep),
      fit_rmst = unname(fit_rmst_val),
      ci_lower = unname(fit_rmst_val - z_scale * se_rmst_val),
      ci_upper = unname(fit_rmst_val + z_scale * se_rmst_val),
    )
}

#' @details The `ertte_coxph` method delegates to
#' `survival::survfit(object, newdata, conf.int = conf_level, se.fit =
#' TRUE)`, the same call [ertte_predict.ertte_coxph()] uses. Because the
#' fitted baseline hazard (and therefore every covariate-adjusted
#' survival curve) is a right-continuous step function, `fit_rmst` is an
#' *exact* finite sum of rectangle areas between consecutive jump times
#' up to `tau` -- not a numerical-quadrature approximation. `se_rmst`
#' comes from a delta method that reuses the same rectangle/tail-weighted
#' construction as `survival:::survmean()`'s classic Greenwood-based RMST
#' variance, but with the variance-increment term replaced by increments
#' of the profile-specific `std.err(t)^2` returned by `survfit()` --
#' `survmean()`'s own Greenwood term is based on population-level risk
#' sets shared across every covariate profile, which understates
#' uncertainty for a profile far from the mean covariate values (see
#' `.ertte_rmst_pfun_delta()`'s source comments for the derivation and
#' the bootstrap cross-check that motivated this).
#'
#' If any value of `tau` exceeds the last observed follow-up time across
#' the whole fitted cohort, `ertte_rmst()` warns: RMST integrates the
#' *entire* curve up to `tau`, so silently assuming survival stays flat
#' beyond the observed range (the same extrapolation convention
#' `ertte_predict.ertte_coxph()` uses for a single time point) has a
#' larger effect on an area than on a point-in-time prediction.
#'
#' `conf_level = 0`/`1`, documented (see `.ertte_check_conf_level()`) as
#' legitimate degenerate endpoints, are supported here directly, since
#' the delta-method interval is built from `z_scale` (`qnorm()`-derived,
#' `0` or `Inf` at these boundaries) rather than `survfit()`'s own
#' `conf.int` machinery -- the latter is only used to request `$surv`/
#' `$std.err`, which don't depend on the requested confidence level (see
#' issue #11); `survfit()` is always called with a fixed, valid
#' placeholder value internally.
#'
#' A zero-row `newdata` returns a zero-row tibble with the expected
#' columns rather than erroring: `survival::survfit()` itself rejects an
#' entirely-missing `newdata` with a cryptic "all rows of newdata have
#' missing values" error (see issue #10).
#'
#' @rdname ertte_rmst
#' @export
#' @examples
#' mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
#' ertte_rmst(mod_cox, ertte_data[1:5, ], tau = c(60, 90))
#'
ertte_rmst.ertte_coxph <- function(object, newdata = NULL, tau, conf_level = .95, ...) {
  .ertte_check_coxph_nevent(object)
  .ertte_check_conf_level(conf_level)
  if (is.null(newdata)) newdata <- object$data
  if (!is.numeric(tau) || length(tau) == 0L || anyNA(tau) || any(tau <= 0)) {
    rlang::abort("`tau` must be a numeric vector of strictly positive values.")
  }
  if (nrow(newdata) == 0L) {
    return(
      newdata |>
        tibble::as_tibble() |>
        dplyr::mutate(
          tau = numeric(0),
          fit_rmst = numeric(0),
          ci_lower = numeric(0),
          ci_upper = numeric(0)
        )
    )
  }
  # `conf.int` here is a fixed, valid placeholder, not `conf_level`
  # itself: only `sf$time`/`sf$surv`/`sf$std.err` are used below (never
  # `sf$lower`/`sf$upper`), and those don't depend on `conf.int` at all
  # -- confirmed empirically. `conf_level` is applied afterwards, via
  # `z_scale`, in `.ertte_rmst_pfun_delta()`'s Wald construction, which
  # handles 0/1 without going through `survfit()`'s own `conf.int`
  # validation (which rejects exactly 0 or 1).
  sf <- survival::survfit(object, newdata = newdata, conf.int = 0.95, se.fit = TRUE)
  max_obs_time <- max(sf$time)
  if (any(tau > max_obs_time)) {
    rlang::warn(paste0(
      "`tau` exceeds the last observed follow-up time (", max_obs_time, ") for ",
      "at least one value. RMST beyond that point assumes survival stays flat ",
      "at its last estimated value, which may understate/overstate the true area."
    ))
  }
  z_scale <- -stats::qnorm((1 - conf_level) / 2)
  start.time <- if (!is.null(sf$t0)) sf$t0 else min(0, sf$time)
  n <- nrow(newdata)
  k <- length(tau)

  out_rows <- vector("list", n * k)
  idx <- 1L
  for (i in seq_len(n)) {
    surv_i <- if (is.matrix(sf$surv)) sf$surv[, i] else sf$surv
    se_i <- if (is.matrix(sf$std.err)) sf$std.err[, i] else sf$std.err
    for (j in seq_len(k)) {
      out_rows[[idx]] <- .ertte_rmst_pfun_delta(sf$time, surv_i, se_i, start.time, tau[j])
      idx <- idx + 1L
    }
  }
  res_mat <- do.call(rbind, out_rows)

  rep_rows <- rep(seq_len(n), each = k)
  tau_rep <- rep(tau, times = n)
  # named distinctly from the `fit_rmst`/`ci_lower`/`ci_upper` columns the
  # `mutate()` below creates -- see the matching comment in
  # `ertte_rmst.ertte_aft()` above (erplots#12).
  fit_rmst_val <- res_mat[, "rmean"]
  se_rmst_val <- res_mat[, "se_rmean"]

  newdata[rep_rows, , drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      tau = unname(tau_rep),
      fit_rmst = unname(fit_rmst_val),
      ci_lower = unname(fit_rmst_val - z_scale * se_rmst_val),
      ci_upper = unname(fit_rmst_val + z_scale * se_rmst_val),
    )
}
