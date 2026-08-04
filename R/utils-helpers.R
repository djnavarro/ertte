
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

.as_ertte_aft <- function(mod, dist) {
  # "ertte_aft" (engine-specific subclass) ahead of "ertte_model" (shared
  # superclass) -- see AGENTS.md "API naming: AFT vs Cox PH" for the
  # dispatch scheme this supports.
  class(mod) <- c("ertte_aft", "ertte_model", class(mod))
  mod$ertte <- list(type = dist) # internal "ertte" list to store ertte-specific info
  mod
}

# All four tested/supported `survreg()` distributions are log-location-scale
# families: log(T) = mu + scale * W, where W has a "base" distribution that
# depends only on `dist`, not on the covariates. This table gives the CDF
# (`pbase`) and quantile function (`qbase`) of that base distribution, used
# by `ertte_predict()` (survival probabilities) and `.ertte_simulate_draws()`
# (inverse-CDF sampling of event times).
.ertte_dist_info <- function(dist) {
  .ertte_check_dist(dist)
  switch(
    dist,
    exponential = ,
    weibull = list(
      pbase = function(z) 1 - exp(-exp(z)),
      qbase = function(p) log(-log(1 - p))
    ),
    lognormal = list(pbase = stats::pnorm, qbase = stats::qnorm),
    loglogistic = list(pbase = stats::plogis, qbase = stats::qlogis)
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
