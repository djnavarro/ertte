
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
# `er_predict.ertte_model()` is currently a minimal placeholder: the
# real landmark-binary / RMST contract described in the package's design
# issue (phase 2 -- scalar E-R views of a TTE endpoint) hasn't been
# designed yet, so this just forwards to `ertte_predict()` with a `time`
# argument threaded through `...`. Revisit once phase 2 is scoped.

er_predict.ertte_model <- function(model, newdata, conf_level = 0.95, ...) {
  dots <- list(...)
  time <- dots$time %||% rlang::abort("er_predict.ertte_model() currently requires a `time` argument (via `...`).")
  ertte_predict(object = model, newdata = newdata, time = time, conf_level = conf_level)
}

er_simulate.ertte_model <- function(model, newdata, nsim = 100, seed = NULL, ...) {
  # `censor_time` isn't part of erplots' fixed `er_simulate(model,
  # newdata, nsim, seed)` contract, but can still be threaded through
  # `...` for callers that want the accurate (rather than default
  # event-rows-uncensored) simulation behaviour -- see
  # `simulate.ertte_model()`'s `censor_time` argument/Details.
  dots <- list(...)
  .ertte_simulate_draws(
    object = model, newdata = newdata, nsim = nsim, seed = seed,
    censor_time = dots$censor_time
  )
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
