
# `conf_level` must be a single number in [0, 1] for the `qnorm()`-based
# z-score used by `ertte_predict()`/`er_summary.ertte_model()` to be
# meaningful. 0 and 1 are legitimate (degenerate) endpoints -- a 0%
# interval collapses to the point estimate (z = 0), and a 100% interval
# is infinitely wide (z = Inf) -- but values outside [0, 1] aren't valid
# probabilities and currently produce a silently reversed or NaN interval
# rather than erroring.
.ertte_check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      conf_level < 0 || conf_level > 1) {
    rlang::abort(paste0(
      "`conf_level` must be a single number between 0 and 1 (inclusive), not ",
      .fmt_bad_value(conf_level), "."
    ))
  }
}

# `term` (as passed to `ertte_add_term()`/`ertte_remove_term()`) must be
# a one-sided formula naming exactly one covariate, e.g. `~ sex` -- same
# validation as erglm's `.erglm_check_term()`.
.ertte_check_term <- function(term) {
  if (is.null(term) || !inherits(term, "formula")) {
    rlang::abort(paste0(
      "`term` must be a one-sided formula naming a single covariate ",
      "(e.g. `~ sex`), not ", .fmt_bad_value(term), "."
    ))
  }
  if (length(term) != 2L) {
    rlang::abort(paste0(
      "`term` must be a one-sided formula (e.g. `~ sex`), not the ",
      "two-sided formula `", deparse(term), "`. ertte_add_term()/",
      "ertte_remove_term() work on plain covariate terms and don't use ",
      "a response."
    ))
  }
  trm_lab <- attr(stats::terms(term), "term.labels")
  if (length(trm_lab) != 1L) {
    rlang::abort(paste0(
      "`term` must name exactly one covariate (e.g. `~ sex`), not ",
      length(trm_lab), ": `", deparse(term), "`."
    ))
  }
}

# `candidates` (as passed to `ertte_scm_forward()`/`ertte_scm_backward()`)
# must be a non-empty character vector where every element names exactly
# one covariate term -- validated up front, same rationale as erglm's
# `.erglm_check_candidates()`.
.ertte_check_candidates <- function(candidates) {
  if (!is.character(candidates) || length(candidates) == 0L || anyNA(candidates)) {
    rlang::abort(paste0(
      "`candidates` must be a non-empty character vector with no missing ",
      "values, not ", .fmt_bad_value(candidates), "."
    ))
  }
  for (cc in candidates) {
    add <- tryCatch(stats::as.formula(paste("~", cc)), error = function(e) NULL)
    if (is.null(add)) {
      rlang::abort(paste0(
        "`candidates` contains an invalid entry: \"", cc, "\" could not ",
        "be parsed as a formula term."
      ))
    }
    trm_lab <- attr(stats::terms(add), "term.labels")
    if (length(trm_lab) != 1L) {
      rlang::abort(paste0(
        "`candidates` contains an invalid entry: \"", cc, "\" names ",
        length(trm_lab), " terms, not exactly one. Each element of ",
        "`candidates` must name a single covariate term (e.g. \"sex\", ",
        "not \"sex + dose\")."
      ))
    }
  }
}

# `nsim` must be a single positive whole number.
.ertte_check_nsim <- function(nsim) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim < 1 || abs(nsim - round(nsim)) > .Machine$double.eps^0.5) {
    rlang::abort(paste0(
      "`nsim` must be a single positive whole number, not ",
      .fmt_bad_value(nsim), "."
    ))
  }
}

# `dist` must name one of the tested/supported AFT distributions.
.ertte_check_dist <- function(dist) {
  supported <- c("exponential", "weibull", "lognormal", "loglogistic")
  if (!is.character(dist) || length(dist) != 1L || is.na(dist) || !dist %in% supported) {
    rlang::abort(paste0(
      "`dist` must be one of \"", paste(supported, collapse = "\", \""),
      "\", not ", .fmt_bad_value(dist), "."
    ))
  }
}

# `time` (as passed to `ertte_predict()` and to the closures returned by
# `ertte_fun()`) must be a numeric vector of strictly positive values --
# `log(time)` (AFT) / the baseline-hazard lookup (Cox) are both undefined
# or meaningless otherwise. Shared by both `ertte_predict.ertte_aft()`/
# `ertte_predict.ertte_coxph()` (which already validated this inline,
# identically, before this helper existed) and both `ertte_fun.ertte_aft()`/
# `ertte_fun.ertte_coxph()`'s returned closures (which previously validated
# nothing, silently returning `NaN`/`1` for non-positive `time` instead of
# erroring -- see issue #7).
.ertte_check_time <- function(time) {
  if (!is.numeric(time) || length(time) == 0L || anyNA(time) || any(time <= 0)) {
    rlang::abort("`time` must be a numeric vector of strictly positive values.")
  }
}

# Extracts the time/event variable names from a two-sided
# `Surv(time, event) ~ ...` formula -- like `.ertte_response_vars()`
# (below), but usable before a model has been fit, since it works
# directly on the formula the caller supplies rather than a fitted
# object's `$terms`. Used by `.ertte_check_response_time()`.
.ertte_surv_vars_from_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    rlang::abort(
      "`formula` must be a two-sided formula with a `survival::Surv()` response."
    )
  }
  lhs <- formula[[2]]
  # allow both `Surv(...)` and `survival::Surv(...)`
  fn_name <- if (is.call(lhs)) deparse(lhs[[1]]) else ""
  if (!grepl("(^|::)Surv$", fn_name)) {
    rlang::abort("The model's response must be a `survival::Surv()` object.")
  }
  args <- as.list(lhs)[-1]
  vars <- vapply(args, deparse, character(1))
  list(time = vars[1], event = vars[2])
}

# `ertte_aft()` inherits `time > 0` validation as a side effect of
# `survival::survreg()`'s own internal check ("Invalid survival times
# for this distribution"), since a log-location-scale AFT model takes
# `log(time)`. `survival::coxph()` has no equivalent check -- a
# negative (or zero) `time` fits silently and produces plausible-
# looking but meaningless downstream predictions (see issue #5).
# Called by `ertte_coxph()` before fitting, so behaviour doesn't depend
# on which underlying engine happens to reject bad input. Only
# non-missing values are checked -- rows with a missing time are
# `coxph()`'s/`survreg()`'s own concern (both drop them silently, the
# base R modelling convention), not this validation's.
.ertte_check_response_time <- function(formula, data) {
  vars <- .ertte_surv_vars_from_formula(formula)
  time_val <- data[[vars$time]]
  if (is.null(time_val)) {
    rlang::abort(paste0(
      "Response variable `", vars$time, "` (from `", deparse(formula[[2]]),
      "`) was not found in `data`."
    ))
  }
  if (!is.numeric(time_val)) {
    rlang::abort(paste0(
      "Response variable `", vars$time, "` must be numeric, not ",
      class(time_val)[1], "."
    ))
  }
  bad <- time_val[!is.na(time_val)] <= 0
  if (any(bad)) {
    rlang::abort(paste0(
      "Response variable `", vars$time, "` must be strictly positive for ",
      "every (non-missing) row -- found ", sum(bad), " non-positive value",
      if (sum(bad) != 1L) "s" else "", "."
    ))
  }
}

# `survival::coxph()` crashes with a cryptic, low-level error when
# fitting on exactly one usable row -- `'x' must be an array of at
# least two dimensions`, from `rowSums(fit$influence[, 1:3])` inside
# `coxph()`'s own post-fit diagnostics, triggered because `influence`
# degrades from a matrix to a plain vector when there's only one
# observation (see issue #4). With zero usable rows, `coxph()` already
# fails with a reasonably clear message of its own ("No (non-missing)
# observations"), so this isn't strictly needed there, but checking
# `< 2` uniformly gives one clear message for both degenerate cases
# rather than depending on which one happens to already be handled
# gracefully upstream.
#
# "Usable rows" means after `coxph()`'s own missing-data handling would
# apply, not simply `nrow(data)` -- a formula with a covariate that is
# `NA` on all-but-one row is just as broken as literally supplying one
# row of data (confirmed empirically). This is approximated by building
# a plain `stats::model.frame()` on `formula`/`data` (with the default
# `na.action`, `na.omit`, matching `coxph()`'s own default) and counting
# its rows -- this only approximates `coxph()`'s real fitting pathway
# (e.g. it ignores `subset`/a caller-supplied `na.action`/`weights` in
# `...`), so it's a heuristic pre-flight check, not a guarantee.
.ertte_check_coxph_data_size <- function(formula, data) {
  mf <- suppressWarnings(
    tryCatch(stats::model.frame(formula, data), error = function(e) NULL)
  )
  n <- if (is.null(mf)) 0L else nrow(mf)
  if (n < 2L) {
    rlang::abort(paste0(
      "ertte_coxph() needs at least 2 usable rows (after removing rows ",
      "with missing values) to fit a Cox model, but only ", n, " row",
      if (n != 1L) "s were" else " was", " available. `survival::coxph()` ",
      "itself fails with a cryptic, low-level error when there's exactly ",
      "one usable row, so this is checked explicitly up front."
    ))
  }
}

# Shared by `.ertte_simulate_draws()`'s per-engine methods: if `seed` is
# `NULL`, pick one and tell the user (since it determines the actual
# simulated values returned), otherwise pass it through unchanged.
.ertte_check_landmark_time <- function(landmark_time) {
  if (!is.numeric(landmark_time) || length(landmark_time) != 1L ||
      is.na(landmark_time) || landmark_time <= 0) {
    rlang::abort(paste0(
      "`landmark_time` must be a single, strictly positive number (a ",
      "fixed time t*), not ", .fmt_bad_value(landmark_time), "."
    ))
  }
}

# `tau`, when threaded through erplots' `er_predict()`/`er_simulate()`
# `...` (see `R/er-methods.R`), must be a single strictly positive
# number: their scalar E-R grammar expects exactly one row (or, for
# simulation, one row per replicate) per `newdata` row, unlike
# `ertte_rmst()`'s own vectorised `tau` argument, which can evaluate
# several horizons at once.
.ertte_check_single_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0) {
    rlang::abort(paste0(
      "`tau` (passed via `...`) must be a single, strictly positive ",
      "number -- erplots' scalar E-R grammar expects one row per ",
      "`newdata` row, so `ertte_rmst()`'s vectorised `tau` argument ",
      "isn't supported through this interface."
    ))
  }
}

.ertte_pick_seed <- function(seed) {
  if (is.null(seed)) {
    seed <- .pick_seed()
    rlang::inform(paste0("Using seed = ", seed, ". Pass `seed = ", seed, "` to reproduce this result."))
  }
  seed
}

# Shared by `.ertte_simulate_draws()`'s per-engine methods: `newdata`
# must contain the original response columns (`time`/`event`, as named
# in `object`'s `Surv()` call) -- used (absent a supplied `censor_time`)
# to cap censored rows' simulated event times at their observed exit
# time, and to identify event rows (left uncensored by default -- see
# `.ertte_apply_admin_censoring()`). Returns the response variable names
# (from `.ertte_response_vars()`) for convenience at the call site.
.ertte_check_newdata_response <- function(object, newdata) {
  vars <- .ertte_response_vars(object)
  if (!all(c(vars$time, vars$event) %in% names(newdata))) {
    rlang::abort(paste0(
      "`newdata` must contain the original response columns `", vars$time,
      "` and `", vars$event, "` (used to cap simulated event times at the ",
      "observed exit time)."
    ))
  }
  vars
}

# Validates `censor_time` (the administrative-censoring-time argument to
# `simulate()`/`.ertte_simulate_draws()`): `NULL`, or a numeric vector of
# strictly positive, non-missing values with length 1 (recycled) or `n`
# (one per row of `newdata`). Returns a length-`n` vector (or `NULL`).
.ertte_check_censor_time <- function(censor_time, n) {
  if (is.null(censor_time)) return(NULL)
  if (!is.numeric(censor_time) || anyNA(censor_time) || any(censor_time <= 0)) {
    rlang::abort(paste0(
      "`censor_time` must be `NULL`, or a numeric vector of strictly ",
      "positive, non-missing values, not ", .fmt_bad_value(censor_time), "."
    ))
  }
  if (!length(censor_time) %in% c(1L, n)) {
    rlang::abort(paste0(
      "`censor_time` must have length 1 (recycled across rows) or ", n,
      " (`nrow(newdata)`), not ", length(censor_time), "."
    ))
  }
  rep_len(censor_time, n)
}

# Shared by `.ertte_simulate_draws()`'s per-engine methods: applies the
# administrative-censoring cap to a vector of raw (uncensored) simulated
# event times, `sim_time_raw`.
#
# If `censor_time` is supplied (already validated/recycled to length
# `n` by `.ertte_check_censor_time()`), it's used as the cap for every
# row uniformly -- the "proper" simulation, whenever a genuine per-row
# (or study-wide constant) administrative follow-up time is known.
#
# Otherwise (the default, `censor_time = NULL`): censored rows
# (`event_obs == 0`) have their simulated draws capped at `obs_time`,
# since that observed exit time *is* the row's true censoring time (the
# mechanism that produced the data really did censor this subject then)
# -- an exact match, not an approximation. Event rows (`event_obs == 1`)
# are left uncensored (`sim_event = 1` always): their observed exit time
# is when the event actually happened, not their administrative
# censoring horizon (which was necessarily later, and is unobserved) --
# capping there would leak the observed event day into the simulated
# draws. This is a documented simplification, still an approximation for
# event rows absent a supplied `censor_time`, but avoids the specific
# bias of the old blanket `pmin(sim_time_raw, obs_time)`.
.ertte_apply_admin_censoring <- function(sim_time_raw, obs_time, event_obs, censor_time) {
  if (!is.null(censor_time)) {
    return(list(
      sim_time = pmin(sim_time_raw, censor_time),
      sim_event = as.numeric(sim_time_raw <= censor_time)
    ))
  }
  is_censored_row <- event_obs == 0
  list(
    sim_time = ifelse(is_censored_row, pmin(sim_time_raw, obs_time), sim_time_raw),
    sim_event = ifelse(is_censored_row, as.numeric(sim_time_raw <= obs_time), 1)
  )
}

.as_ertte_aft <- function(mod, dist) {
  # "ertte_aft" (engine-specific subclass) ahead of "ertte_model" (shared
  # superclass) -- see AGENTS.md "API naming: AFT vs Cox PH" for the
  # dispatch scheme this supports.
  class(mod) <- c("ertte_aft", "ertte_model", class(mod))
  mod$ertte <- list(type = dist) # internal "ertte" list to store ertte-specific info
  mod
}

.as_ertte_coxph <- function(mod) {
  # "ertte_coxph" (engine-specific subclass) ahead of "ertte_model"
  # (shared superclass) -- see AGENTS.md "API naming: AFT vs Cox PH" for
  # the dispatch scheme this supports. `type` is fixed at `"coxph"`
  # (there's no `dist` equivalent to select, unlike `.as_ertte_aft()`),
  # kept mainly so code that reads `object$ertte$type` (e.g.
  # `er_summary.ertte_model()`) has something sensible to report.
  class(mod) <- c("ertte_coxph", "ertte_model", class(mod))
  mod$ertte <- list(type = "coxph")
  mod
}

# All four tested/supported `survreg()` distributions are log-location-scale
# families: log(T) = mu + scale * W, where W has a "base" distribution that
# depends only on `dist`, not on the covariates. This table gives the CDF
# (`pbase`), quantile function (`qbase`), and density (`dbase`) of that base
# distribution. `pbase`/`qbase` are used by `ertte_predict()` (survival
# probabilities) and `.ertte_simulate_draws()` (inverse-CDF sampling of
# event times); `dbase` is used by `ertte_rmst.ertte_aft()` (its delta-method
# gradient, `d/dmu S(t|x) = dbase(z) / scale`, needs the base density).
.ertte_dist_info <- function(dist) {
  .ertte_check_dist(dist)
  switch(
    dist,
    exponential = ,
    weibull = list(
      pbase = function(z) 1 - exp(-exp(z)),
      qbase = function(p) log(-log(1 - p)),
      dbase = function(z) exp(z - exp(z))
    ),
    lognormal = list(pbase = stats::pnorm, qbase = stats::qnorm, dbase = stats::dnorm),
    loglogistic = list(pbase = stats::plogis, qbase = stats::qlogis, dbase = stats::dlogis)
  )
}

# Extracts the names of the time/event variables from the two-sided
# `Surv(time, event) ~ ...` formula a fitted `ertte_model` was built from --
# used by `.ertte_simulate_draws()` to find the observed censoring/follow-up
# time to cap simulated event times at.
.ertte_response_vars <- function(object) {
  lhs <- object$terms[[2]]
  # allow both `Surv(...)` and `survival::Surv(...)`
  fn_name <- if (is.call(lhs)) deparse(lhs[[1]]) else ""
  if (!grepl("(^|::)Surv$", fn_name)) {
    rlang::abort("The model's response must be a `survival::Surv()` object.")
  }
  args <- as.list(lhs)[-1]
  vars <- vapply(args, deparse, character(1))
  list(time = vars[1], event = vars[2])
}

.ertte_term_in_model <- function(mod, term) {
  trm_mod <- stats::terms(mod)
  trm_tst <- stats::terms(term)
  trm_mod_lab <- attr(trm_mod, "term.labels")
  trm_tst_lab <- attr(trm_tst, "term.labels")
  ind <- which(trm_mod_lab == trm_tst_lab)
  return(length(ind) != 0)
}
