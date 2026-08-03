
#' Stepwise covariate modelling for exposure-response TTE models
#'
#' @param mod An ertte model object
#' @param candidates Character vector with list of candidate terms
#' @param threshold Threshold to test against
#' @param seed Optional seed to control order of term tests
#'
#' @returns For `ertte_scm_forward()` and `ertte_scm_backward()`, the
#' updated ertte model is returned, with the SCM history log updated
#' internally. For `ertte_scm_history()`, a data frame is returned
#' containing the SCM history log
#'
#' @details Terms are compared with a likelihood-ratio Chi-squared test
#' (`stats::anova()` on nested `survreg` fits) -- unlike the companion
#' `erglm` package's SCM, there's no family-dependent choice of test
#' here, since a `survreg` model's likelihood ratio test doesn't vary by
#' distribution.
#'
#' `seed` exists as a safety measure against run-to-run variation in the
#' order candidate terms are tested within a step (`sample()`, shuffled
#' before testing one at a time). Model fitting itself
#' (`survival::survreg()`) is deterministic given a starting formula, so
#' `seed` only matters in the (essentially measure-zero) case of an
#' exact p-value tie between competing candidates within a step -- see
#' the companion `erglm` package's equivalent documentation for the full
#' rationale, which applies unchanged here.
#'
#' If a candidate term is aliased (perfectly collinear) with a term
#' already in the model, `stats::anova()` reports an `NA` p-value for
#' it. That candidate is skipped for the step (with a warning) rather
#' than being selected or crashing the search.
#'
#' `candidates` is validated up front: every element must be parseable
#' as a formula and name exactly one covariate term (e.g. `"sex"`, not
#' `"sex + dose"` or `"not a formula"`).
#'
#' @name ertte_scm
#' @examples
#' mod0 <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
#' mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"))
#' ertte_scm_history(mod1)
#'
#' mod2 <- ertte_model(survival::Surv(time, event) ~ aucss + sex + dose, ertte_data)
#' mod3 <- ertte_scm_backward(mod2, candidates = c("sex", "dose"))
#' ertte_scm_history(mod3)
NULL

#' @rdname ertte_scm
#' @export
ertte_scm_forward <- function(mod, candidates, threshold = 0.01, seed = NULL) {
  .ertte_check_candidates(candidates)
  if (is.null(seed)) {
    seed <- .pick_seed()
  }
  withr::with_seed(
    seed = seed,
    code = {
      mod_out <- .ertte_scm_forward(mod = mod, candidates = candidates, threshold = threshold)
    }
  )
  return(mod_out)
}

.ertte_scm_forward <- function(mod, candidates, threshold) {
  history <- ertte_scm_history(mod)
  last_iter <- max(history$iteration)
  while (TRUE) {
    mod_new <- .ertte_once_forward(mod, candidates, threshold)
    history_new <- ertte_scm_history(mod_new)
    this_iter <- max(history_new$iteration)
    if (this_iter == last_iter) return(mod)
    history <- history_new
    last_iter <- this_iter
    mod <- mod_new
    updates <- history |>
      dplyr::filter(iteration == last_iter) |>
      dplyr::pull(model_updated)
    if (all(updates == 0L)) return(mod)
  }
}

#' @rdname ertte_scm
#' @export
ertte_scm_backward <- function(mod, candidates, threshold = 0.001, seed = NULL) {
  .ertte_check_candidates(candidates)
  if (is.null(seed)) {
    seed <- .pick_seed()
  }
  withr::with_seed(
    seed = seed,
    code = {
      mod_out <- .ertte_scm_backward(mod = mod, candidates = candidates, threshold = threshold)
    }
  )
  return(mod_out)
}

.ertte_scm_backward <- function(mod, candidates, threshold) {
  history <- ertte_scm_history(mod)
  last_iter <- max(history$iteration)
  while (TRUE) {
    mod_new <- .ertte_once_backward(mod, candidates, threshold)
    history_new <- ertte_scm_history(mod_new)
    this_iter <- max(history_new$iteration)
    if (this_iter == last_iter) return(mod)
    history <- history_new
    last_iter <- this_iter
    mod <- mod_new
    updates <- history |>
      dplyr::filter(iteration == last_iter) |>
      dplyr::pull(model_updated)
    if (all(updates == 0L)) return(mod)
  }
}

#' @rdname ertte_scm
#' @export
ertte_scm_history <- function(mod) {
  history <- mod$ertte$history
  if (!is.null(history)) return(history)
  history_row <- tibble::tibble(
    iteration = 0L,
    attempt = 0L,
    step = "base model",
    action = NA_character_,
    term_tested = NA_character_,
    model_tested = deparse(stats::formula(mod)),
    model_converged = isTRUE(mod$converged) || is.null(mod$converged),
    term_p_value = NA_real_,
    model_aic = stats::AIC(mod),
    model_bic = stats::BIC(mod),
    model_updated = NA
  )
  return(history_row)
}

.ertte_once_forward <- function(mod, candidates, threshold) {
  candidates <- sample(candidates)
  history <- ertte_scm_history(mod)
  iter <- max(history$iteration) + 1L
  attm <- max(history$attempt)
  lowest_p <- threshold
  update_ind <- NA_integer_
  best_mod <- mod
  for (cc in candidates) {
    add <- stats::as.formula(paste("~", cc))
    attm <- attm + 1L
    if (!.ertte_term_in_model(mod, add)) {
      mod_new <- ertte_add_term(mod, add, quiet = TRUE)
      p_val <- .ertte_anova_p(mod, mod_new)
      history_row <- tibble::tibble(
        iteration = iter,
        attempt = attm,
        step = "forward",
        action = "add",
        term_tested = deparse(add),
        model_tested = deparse(stats::formula(mod_new)),
        model_converged = isTRUE(mod_new$converged) || is.null(mod_new$converged),
        term_p_value = p_val,
        model_aic = stats::AIC(mod_new),
        model_bic = stats::BIC(mod_new),
        model_updated = NA
      )
      history <- tibble::add_row(history, history_row)
      if (is.na(p_val)) {
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(add), "` in forward step ",
          iter, ": comparison p-value is NA (often caused by a candidate ",
          "that's aliased/collinear with a term already in the model)."
        ))
      } else if (p_val < lowest_p) {
        update_ind <- attm
        lowest_p <- p_val
        best_mod <- mod_new
      }
    }
  }
  history <- history |>
    dplyr::mutate(
      model_updated = dplyr::case_when(
        iteration != iter ~ model_updated,
        attempt == update_ind ~ 1L,
        TRUE ~ 0L
      )
    )
  best_mod$ertte$history <- history
  return(best_mod)
}

.ertte_once_backward <- function(mod, candidates, threshold) {
  trm_mod <- stats::terms(mod)
  trm_lab <- attr(trm_mod, "term.labels")
  candidates <- intersect(trm_lab, candidates)
  if (length(candidates) == 0L) return(mod)
  candidates <- sample(candidates)
  history <- ertte_scm_history(mod)
  iter <- max(history$iteration) + 1L
  attm <- max(history$attempt)
  highest_p <- threshold
  update_ind <- NA_integer_
  best_mod <- mod
  for (cc in candidates) {
    del <- stats::as.formula(paste("~", cc))
    attm <- attm + 1L
    if (.ertte_term_in_model(mod, del)) {
      mod_new <- ertte_remove_term(mod, del, quiet = TRUE)
      p_val <- .ertte_anova_p(mod, mod_new)
      history_row <- tibble::tibble(
        iteration = iter,
        attempt = attm,
        step = "backward",
        action = "remove",
        term_tested = deparse(del),
        model_tested = deparse(stats::formula(mod_new)),
        model_converged = isTRUE(mod_new$converged) || is.null(mod_new$converged),
        term_p_value = p_val,
        model_aic = stats::AIC(mod_new),
        model_bic = stats::BIC(mod_new),
        model_updated = NA
      )
      history <- tibble::add_row(history, history_row)
      if (is.na(p_val)) {
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(del), "` in backward step ",
          iter, ": comparison p-value is NA (often caused by a candidate ",
          "that's aliased/collinear with another term in the model)."
        ))
      } else if (p_val > highest_p) {
        update_ind <- attm
        highest_p <- p_val
        best_mod <- mod_new
      }
    }
  }
  history <- history |>
    dplyr::mutate(
      model_updated = dplyr::case_when(
        iteration != iter ~ model_updated,
        attempt == update_ind ~ 1L,
        TRUE ~ 0L
      )
    )
  best_mod$ertte$history <- history
  return(best_mod)
}

.ertte_anova_p <- function(mod1, mod2) {
  smm <- stats::anova(mod1, mod2)
  p_col <- grep("^Pr\\(", colnames(smm))[1]
  return(smm[[p_col]][2])
}

#' Add or remove a covariate term from an exposure-response TTE model
#'
#' Add or remove a single covariate term from an existing ertte model,
#' returning a new fitted model object.
#'
#' @param mod An ertte model object, as returned by [ertte_model()]
#' @param term A one-sided formula naming the term to add/remove, e.g.
#' `~ sex`
#' @param quiet If `TRUE`, suppress the warning issued when the term
#' can't be added/removed (because it's already in the model / isn't in
#' the model, respectively)
#'
#' @details These functions are not typically called directly; they
#' underpin [ertte_scm_forward()] and [ertte_scm_backward()]. Named and
#' shaped to match the companion `erglm` package's
#' `erglm_add_term()`/`erglm_remove_term()`: covariates enter as plain
#' additive terms on the AFT location scale (linear for continuous
#' covariates, factor levels for categorical ones) -- the richer
#' "continuous covariate as power function" parameterisation described
#' in the package's design issue is not yet implemented (see AGENTS.md).
#'
#' @returns An ertte model object. If the term can't be added/removed
#' (see `quiet`), the original `mod` is returned unchanged.
#'
#' @name ertte_term
#' @examples
#' mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
#' mod2 <- ertte_add_term(mod, ~ sex)
#' mod3 <- ertte_remove_term(mod2, ~ sex)
NULL

#' @rdname ertte_term
#' @export
ertte_add_term <- function(mod, term, quiet = FALSE) {
  .ertte_check_term(term)
  trm_mod <- stats::terms(mod)
  trm_add <- stats::terms(term)
  trm_mod_lab <- attr(trm_mod, "term.labels")
  trm_add_lab <- attr(trm_add, "term.labels")
  ind <- which(trm_mod_lab == trm_add_lab)
  if (length(ind) != 0L) {
    if (!quiet) rlang::warn("cannot add a term that already exists in the model")
    return(mod)
  }
  trm_add_var <- all.vars(attr(trm_add, "variables"))
  dat <- mod$data
  vars_ok <- trm_add_var %in% names(dat)
  if (!all(vars_ok)) {
    if (!quiet) rlang::warn("cannot add a term that uses variables not in the data")
    return(mod)
  }
  fml <- stats::as.formula(
    paste(deparse(stats::formula(mod)), deparse(term[[2]]), sep = " + ")
  )
  ertte_model(formula = fml, data = dat, dist = mod$ertte$type)
}

#' @rdname ertte_term
#' @export
ertte_remove_term <- function(mod, term, quiet = FALSE) {
  .ertte_check_term(term)
  trm_mod <- stats::terms(mod)
  trm_del <- stats::terms(term)
  trm_mod_lab <- attr(trm_mod, "term.labels")
  trm_del_lab <- attr(trm_del, "term.labels")
  ind <- which(trm_mod_lab == trm_del_lab)
  if (length(ind) == 0L) {
    if (!quiet) rlang::warn("cannot remove a term that does not exist in the model")
    return(mod)
  }
  dat <- mod$data
  trm_new <- stats::drop.terms(trm_mod, ind, keep.response = TRUE)
  ertte_model(formula = trm_new, data = dat, dist = mod$ertte$type)
}
