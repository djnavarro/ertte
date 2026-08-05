
#' Select an AFT distribution by AIC
#'
#' Fits each candidate AFT distribution to the same formula/data and
#' selects the best fit by AIC.
#'
#' @param formula Model formula, as for [ertte_aft()]
#' @param data Data set
#' @param candidates Character vector of candidate `dist` values to try.
#' Defaults to all four tested/supported distributions.
#' @returns A list with two elements: `comparison` (a tibble with one row
#' per candidate: `dist`, `logLik`, `aic`, `bic`, `converged`, sorted by
#' AIC) and `model` (the best-fitting `ertte_aft` model, i.e. the one
#' with lowest AIC).
#'
#' @details Ties (to floating point) are broken by the order `candidates`
#' is given in, i.e. the first-listed of the tied candidates is
#' returned as `model`.
#'
#' `candidates` must be a non-empty character vector with no missing
#' values; each element must also separately name one of the
#' tested/supported distributions (see [ertte_aft()]'s `dist`
#' argument). `candidates = character(0)` errors rather than silently
#' fitting nothing and returning a degenerate `list(comparison = <0-row
#' tibble>, model = NULL)`.
#'
#' @export
#' @examples
#' # (uses ertte_aft() internally for each candidate distribution)
#' cmp <- ertte_aft_select_distribution(Surv(time, event) ~ aucss, ertte_data)
#' cmp$comparison
#' cmp$model
#'
ertte_aft_select_distribution <- function(formula, data, candidates = c("exponential", "weibull", "lognormal", "loglogistic")) {
  .ertte_check_dist_candidates(candidates)
  for (cc in candidates) .ertte_check_dist(cc)
  fits <- lapply(candidates, function(dd) ertte_aft(formula, data, dist = dd))
  names(fits) <- candidates
  comparison <- tibble::tibble(
    dist = candidates,
    logLik = vapply(fits, function(mm) as.numeric(stats::logLik(mm)), numeric(1)),
    aic = vapply(fits, stats::AIC, numeric(1)),
    bic = vapply(fits, stats::BIC, numeric(1)),
    converged = vapply(fits, function(mm) isTRUE(mm$converged) || is.null(mm$converged), logical(1)),
  ) |>
    dplyr::arrange(aic)
  best_dist <- comparison$dist[1]
  list(comparison = comparison, model = fits[[best_dist]])
}
