
#' Power-function covariate transform for exposure-response TTE models
#'
#' A formula helper used to enter a continuous covariate as a power 
#' transformed covariate rather than a plain linear term.
#'
#' @param x A numeric covariate. Every non-missing value must be strictly
#' positive.
#' @param ref A single strictly positive reference value. Defaults to
#' the median value of `x`.
#'
#' @returns A numeric vector equal to `log(x / ref)`, classed as an
#' `"ertte_power"` object, with the value of `ref` stored as an attribute.
#'
#' @details 
#' Both [ertte_aft()] and [ertte_coxph()] are linear in their covariates on 
#' the model's natural scale (i.e., log-time for AFT models and log-hazard-ratio
#' for Cox models). A power-function covariate effect is defined as follows for
#' an AFT model:
#' 
#' \deqn{t = t_\mbox{ref} \times (x / x_\mbox{ref})^\theta}
#' 
#' The equivalent for Cox models is as follows:
#' 
#' \deqn{h(t | x) = h_0(t) \times (x / x_\mbox{ref})^\theta}
#' 
#' Once log-transformed, these become linear terms in `log(x / ref)`.
#'
#' So `ertte_power(x)` reduces the power-function parameterisation to an
#' ordinary covariate column: the fitted `survreg()`/`coxph()` coefficient
#' on `ertte_power(x)` *is* the power exponent `theta` directly, and its
#' ordinary Wald confidence interval (from `confint()`/`summary()`) is
#' exactly the confidence interval on `theta`.
#'
#' `ertte_power()` requires every non-missing value of `x` to be strictly
#' positive, so cannot be used for covariates with a placebo/zero-dose
#' group (e.g. `dose`, `aucss`, `cmaxss` in `ertte_data`). It is not 
#' generally appropriate for the exposure itself, only for (some) covariate
#' terms in the TTE model.
#'
#' Nothing prevents combining a plain linear term (`age`) and a power term
#' (`ertte_power(age)`) for the same underlying variable in the same model
#' or SCM candidate set -- term handling works on formula term-labels, not
#' variable semantics, so this is left to the user's judgement.
#'
#' `ref` is fixed at fitting time from the data `ertte_power()` is
#' evaluated on, and reused (not recomputed) when the fitted model is used
#' to predict or simulate on new data, via an internal
#' `stats::makepredictcall()` method -- the same mechanism `stats::poly()`
#' uses for this purpose.
#' 
#' @export
#' @examples
#' mod <- ertte_aft(Surv(time, event) ~ aucss + ertte_power(age), ertte_data)
#' summary(mod)
#'
#' # reference value used for the power transform
#' attr(ertte_power(ertte_data$age), "ref")
#' 
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

#' @exportS3Method stats::makepredictcall
makepredictcall.ertte_power <- function(var, call) {
  if (as.character(call)[1L] == "ertte_power") {
    call$ref <- attr(var, "ref")
  }
  call
}
