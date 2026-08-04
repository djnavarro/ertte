
#' Power-function covariate transform for exposure-response TTE models
#'
#' A `poly()`/`splines::ns()`-style formula helper for entering a
#' continuous covariate as a power function, e.g. `ertte_power(age)` in
#' place of a plain `age` term.
#'
#' @param x A numeric covariate. Every non-missing value must be strictly
#' positive (`log(x / ref)` is undefined otherwise).
#' @param ref A single strictly positive reference value. Defaults to
#' `median(x, na.rm = TRUE)` -- the usual pop-PK/NONMEM convention of
#' referencing the covariate's typical (median) value in the fitting data.
#'
#' @returns A numeric vector equal to `log(x / ref)`, classed
#' `"ertte_power"`, with `ref` stored as an attribute.
#'
#' @details Both [ertte_aft()] (a log-location-scale AFT model: `log(T) =
#' mu + scale * W`) and `ertte_coxph()` (a Cox PH model: `h(t | x) = h0(t)
#' * exp(x'beta)`) are linear in their covariates on the model's natural
#' (log-time or log-hazard-ratio) scale. A power-function covariate effect
#' -- `T = T_ref * (x / ref)^theta` on the AFT time scale, or `h(t | x) =
#' h0(t) * (x / ref)^theta` on the Cox hazard scale -- is, after taking
#' logs, exactly a linear term in `log(x / ref)`:
#'
#' ```
#' log(T)    = ... + theta * log(x / ref) + scale * W
#' log(h/h0) = ... + theta * log(x / ref)
#' ```
#'
#' So `ertte_power(x)` reduces the power-function parameterisation to an
#' ordinary covariate column: the fitted `survreg()`/`coxph()` coefficient
#' on `ertte_power(x)` *is* the power exponent `theta` directly, and its
#' ordinary Wald confidence interval (from `confint()`/`summary()`) is
#' exactly the confidence interval on `theta` -- no delta method or
#' profile-likelihood machinery is needed, unlike covariate power
#' functions on genuinely nonlinear structural parameters (e.g. the
#' companion `emaxnls` package's Emax/EC50 parameters, where this
#' reduction doesn't apply).
#'
#' `ref` is fixed at fitting time from the data `ertte_power()` is
#' evaluated on, and reused (not recomputed) when the fitted model is used
#' to predict on new data -- via a `makepredictcall.ertte_power()` method,
#' the same mechanism `stats::poly()`/`splines::ns()` use for this purpose.
#'
#' `ertte_power()` requires every non-missing value of `x` to be strictly
#' positive, which rules it out for covariates with a placebo/zero-dose
#' group (e.g. `dose`, `aucss`, `cmaxss` in `ertte_data`). This is by
#' design: the power-function parameterisation described in the package's
#' design issue is aimed at the *covariate model* (e.g. age, weight), not
#' the primary exposure metric, which enters the model directly (see
#' [ertte_aft()]).
#'
#' Nothing in [ertte_add_term()]/[ertte_scm_forward()]/
#' [ertte_scm_backward()] prevents combining a plain linear term (`age`)
#' and a power term (`ertte_power(age)`) for the same underlying variable
#' -- term handling throughout ertte works on formula term-labels, not
#' variable semantics, so this is left to the user's judgement.
#'
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss + ertte_power(age), ertte_data)
#' summary(mod)
#'
#' # reference value used for the power transform
#' attr(ertte_power(ertte_data$age), "ref")
ertte_power <- function(x, ref = NULL) {
  if (!is.numeric(x)) {
    rlang::abort(paste0("`x` must be numeric, not ", .fmt_bad_value(x), "."))
  }
  if (any(!is.na(x) & x <= 0)) {
    rlang::abort(paste0(
      "ertte_power() requires every non-missing value of `x` to be ",
      "strictly positive (`log(x / ref)` is undefined otherwise)."
    ))
  }
  if (is.null(ref)) {
    ref <- stats::median(x, na.rm = TRUE)
  } else if (!is.numeric(ref) || length(ref) != 1L || is.na(ref) || ref <= 0) {
    rlang::abort(paste0(
      "`ref` must be a single strictly positive number, not ",
      .fmt_bad_value(ref), "."
    ))
  }
  out <- log(x / ref)
  attr(out, "ref") <- ref
  class(out) <- c("ertte_power", class(out))
  out
}

#' @details `makepredictcall.ertte_power()` is a `stats::makepredictcall()`
#' method, not typically called directly. It ensures that when a fitted
#' model containing an `ertte_power()` term is used to predict/simulate on
#' new data (via `stats::model.matrix()`/`stats::model.frame()` on the
#' model's `terms()`), the *original* fitting-time `ref` is reused rather
#' than a new one recomputed from whatever data is supplied -- the same
#' mechanism `stats::poly()`/`splines::ns()` use.
#'
#' @param var The evaluated variable (here, the `ertte_power()`-transformed
#' vector from the original model fit).
#' @param call The unevaluated call to be reconstructed for new data.
#'
#' @rdname ertte_power
#' @importFrom stats makepredictcall
#' @export
makepredictcall.ertte_power <- function(var, call) {
  if (as.character(call)[1L] == "ertte_power") {
    call$ref <- attr(var, "ref")
  }
  call
}
