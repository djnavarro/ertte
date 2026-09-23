
#' Stepwise covariate modelling for exposure-response TTE models
#'
#' Iteratively adds (`ertte_scm_forward()`) or removes
#' (`ertte_scm_backward()`) covariate terms from an ertte model one at a
#' time, selecting each step's term by a likelihood-ratio p-value, AIC, or
#' BIC, and logging the search in a history object retrievable via
#' `ertte_scm_history()`.
#'
#' @param mod An ertte model object
#' @param candidates Character vector with list of candidate terms
#' @param threshold Threshold to test against. Used only when
#' `criterion = "p-value"` (the default); ignored otherwise. Defaults to
#' `0.01` for `ertte_scm_forward()` and `0.001` for `ertte_scm_backward()`.
#' @param criterion Model selection criterion. One of `"p-value"`
#' (default), `"aic"`, or `"bic"`.
#' @param seed Optional seed to control the order candidate terms are
#' tested within a step. If `NULL` (the default), one is chosen
#' automatically (see "Candidate test order and `seed`" below for when
#' this actually matters).
#'
#' @return 
#' For `ertte_scm_forward()` and `ertte_scm_backward()`, the
#' updated ertte model is returned, with the SCM history log updated
#' internally. For `ertte_scm_history()`, a data frame is returned
#' containing the SCM history log.
#'
#' @details
#' Three model selection criteria are available via the `criterion`
#' argument:
#'
#' - `"p-value"` (default): Models are compared with a likelihood-ratio 
#'   Chi-squared test. A term is added if its likelihood-ratio
#'   p-value falls below `threshold` (forward) or removed if its p-value
#'   exceeds `threshold` (backward). When multiple candidates satisfy the
#'   threshold within a step, the one with the most extreme p-value is
#'   chosen.
#' - `"aic"`: A term is added (forward) or removed (backward) if doing so
#'   strictly decreases AIC relative to the current model. When multiple
#'   candidates improve AIC, the one yielding the lowest AIC is chosen.
#' - `"bic"`: Same as `"aic"`, but using BIC as the criterion.
#'
#' When `criterion` is `"aic"` or `"bic"`, the `threshold` argument has no
#' effect and is ignored, and `term_p_value` is left `NA` in the history
#' for every candidate tested that step (the likelihood-ratio test isn't
#' computed, since it plays no role in selection). The `model_aic` and 
#' `model_bic` columns are always recorded regardless of which criterion 
#' drove selection, and the history's `criterion` column records which one 
#' was used for each forward/backward step.
#'
#' @section Handling problem candidates:
#' A candidate that can't be fit for a given step (e.g. it's aliased with
#' a term already in the model, references a variable missing from the
#' data, or fails to converge) is skipped for that step with a warning,
#' rather than aborting the whole search or being silently selected.
#'
#' @section Candidate test order and `seed`:
#' `seed` exists as a safety measure against run-to-run variation in the
#' order candidate terms are tested within a step (`sample()`, shuffled
#' before testing one at a time). Model fitting itself
#' (`survival::survreg()`) is deterministic given a starting formula, so
#' `seed` only matters in the (essentially measure-zero) case of an
#' exact p-value tie between competing candidates within a step -- see
#' the companion `erglm` package's equivalent documentation for the full
#' rationale, which applies unchanged here.
#'
#' @name ertte_scm
#' @examples
#' # Forward addition by p-value selection
#' mod0 <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"))
#' ertte_scm_history(mod1)
#'
#' # Backward elimination by p-value selection
#' mod2 <- ertte_aft(Surv(time, event) ~ aucss + sex + dose, ertte_data)
#' mod3 <- ertte_scm_backward(mod2, candidates = c("sex", "dose"))
#' ertte_scm_history(mod3)
#'
#' # AIC-based forward addition/backward elimination instead of p-value
#' mod4 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"), criterion = "aic")
#' mod5 <- ertte_scm_backward(mod4, candidates = c("sex", "dose"), criterion = "bic")
#' ertte_scm_history(mod5)
NULL

#' @rdname ertte_scm
#' @export
ertte_scm_forward <- function(mod, candidates, threshold = 0.01, criterion = "p-value", seed = NULL) {
  .ertte_check_candidates(candidates)
  .ertte_check_criterion(criterion)
  if (is.null(seed)) {
    seed <- .pick_seed()
  }
  withr::with_seed(
    seed = seed,
    code = {
      mod_out <- .ertte_scm_forward(mod = mod, candidates = candidates, threshold = threshold, criterion = criterion)
    }
  )
  return(mod_out)
}

.ertte_scm_forward <- function(mod, candidates, threshold, criterion = "p-value") {
  history <- ertte_scm_history(mod)
  last_iter <- max(history$iteration)
  while (TRUE) {
    mod_new <- .ertte_once_forward(mod, candidates, threshold, criterion)
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
ertte_scm_backward <- function(mod, candidates, threshold = 0.001, criterion = "p-value", seed = NULL) {
  .ertte_check_candidates(candidates)
  .ertte_check_criterion(criterion)
  if (is.null(seed)) {
    seed <- .pick_seed()
  }
  withr::with_seed(
    seed = seed,
    code = {
      mod_out <- .ertte_scm_backward(mod = mod, candidates = candidates, threshold = threshold, criterion = criterion)
    }
  )
  return(mod_out)
}

.ertte_scm_backward <- function(mod, candidates, threshold, criterion = "p-value") {
  history <- ertte_scm_history(mod)
  last_iter <- max(history$iteration)
  while (TRUE) {
    mod_new <- .ertte_once_backward(mod, candidates, threshold, criterion)
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
    criterion = NA_character_,
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

# `use_ic`/`ic_fn` implement the "aic"/"bic" branch of `criterion`,
# mirroring the companion `emaxnls` package's development version: a
# candidate is compared against `best_metric` (the *current* model's IC,
# updated as better candidates are found within the step) rather than
# against `threshold`, which only applies to `criterion = "p-value"`. The
# likelihood-ratio p-value isn't computed at all in IC mode -- it plays
# no role in selection there -- so `term_p_value` is left `NA` in the
# history for every row tested under "aic"/"bic".
.ertte_once_forward <- function(mod, candidates, threshold, criterion = "p-value") {
  candidates <- sample(candidates)
  history <- ertte_scm_history(mod)
  iter <- max(history$iteration) + 1L
  attm <- max(history$attempt)
  use_ic <- criterion %in% c("aic", "bic")
  ic_fn <- if (criterion == "bic") stats::BIC else stats::AIC
  best_metric <- if (use_ic) as.numeric(ic_fn(mod)) else threshold
  update_ind <- NA_integer_
  best_mod <- mod
  for (cc in candidates) {
    add <- stats::as.formula(paste("~", cc))
    attm <- attm + 1L
    if (!.ertte_term_in_model(mod, add)) {
      mod_new <- tryCatch(
        ertte_add_term(mod, add, quiet = TRUE),
        error = function(e) e
      )
      if (inherits(mod_new, "error")) {
        history_row <- tibble::tibble(
          iteration = iter,
          attempt = attm,
          step = "forward",
          criterion = criterion,
          action = "add",
          term_tested = deparse(add),
          model_tested = NA_character_,
          model_converged = FALSE,
          term_p_value = NA_real_,
          model_aic = NA_real_,
          model_bic = NA_real_,
          model_updated = NA
        )
        history <- tibble::add_row(history, history_row)
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(add), "` in forward step ",
          iter, ": refitting the model with this term failed with error: ",
          conditionMessage(mod_new)
        ))
        next
      }
      if (identical(stats::formula(mod_new), stats::formula(mod))) {
        history_row <- tibble::tibble(
          iteration = iter,
          attempt = attm,
          step = "forward",
          criterion = criterion,
          action = "add",
          term_tested = deparse(add),
          model_tested = deparse(stats::formula(mod_new)),
          model_converged = isTRUE(mod_new$converged) || is.null(mod_new$converged),
          term_p_value = NA_real_,
          model_aic = stats::AIC(mod_new),
          model_bic = stats::BIC(mod_new),
          model_updated = NA
        )
        history <- tibble::add_row(history, history_row)
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(add), "` in forward step ",
          iter, ": the term could not be added -- it may reference a ",
          "variable not present in the fitting data. See ertte_add_term() ",
          "for details."
        ))
        next
      }
      p_val <- if (use_ic) NA_real_ else .ertte_anova_p(mod, mod_new)
      history_row <- tibble::tibble(
        iteration = iter,
        attempt = attm,
        step = "forward",
        criterion = criterion,
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
      if (use_ic) {
        candidate_ic <- as.numeric(ic_fn(mod_new))
        if (candidate_ic < best_metric) {
          update_ind <- attm
          best_metric <- candidate_ic
          best_mod <- mod_new
        }
      } else if (is.na(p_val)) {
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(add), "` in forward step ",
          iter, ": comparison p-value is NA (often caused by a candidate ",
          "that's aliased/collinear with a term already in the model)."
        ))
      } else if (p_val < best_metric) {
        update_ind <- attm
        best_metric <- p_val
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

.ertte_once_backward <- function(mod, candidates, threshold, criterion = "p-value") {
  trm_mod <- stats::terms(mod)
  trm_lab <- attr(trm_mod, "term.labels")
  candidates <- intersect(trm_lab, candidates)
  if (length(candidates) == 0L) return(mod)
  candidates <- sample(candidates)
  history <- ertte_scm_history(mod)
  iter <- max(history$iteration) + 1L
  attm <- max(history$attempt)
  use_ic <- criterion %in% c("aic", "bic")
  ic_fn <- if (criterion == "bic") stats::BIC else stats::AIC
  best_metric <- if (use_ic) as.numeric(ic_fn(mod)) else threshold
  update_ind <- NA_integer_
  best_mod <- mod
  for (cc in candidates) {
    del <- stats::as.formula(paste("~", cc))
    attm <- attm + 1L
    if (.ertte_term_in_model(mod, del)) {
      mod_new <- tryCatch(
        ertte_remove_term(mod, del, quiet = TRUE),
        error = function(e) e
      )
      if (inherits(mod_new, "error")) {
        history_row <- tibble::tibble(
          iteration = iter,
          attempt = attm,
          step = "backward",
          criterion = criterion,
          action = "remove",
          term_tested = deparse(del),
          model_tested = NA_character_,
          model_converged = FALSE,
          term_p_value = NA_real_,
          model_aic = NA_real_,
          model_bic = NA_real_,
          model_updated = NA
        )
        history <- tibble::add_row(history, history_row)
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(del), "` in backward step ",
          iter, ": refitting the model without this term failed with error: ",
          conditionMessage(mod_new)
        ))
        next
      }
      if (identical(stats::formula(mod_new), stats::formula(mod))) {
        history_row <- tibble::tibble(
          iteration = iter,
          attempt = attm,
          step = "backward",
          criterion = criterion,
          action = "remove",
          term_tested = deparse(del),
          model_tested = deparse(stats::formula(mod_new)),
          model_converged = isTRUE(mod_new$converged) || is.null(mod_new$converged),
          term_p_value = NA_real_,
          model_aic = stats::AIC(mod_new),
          model_bic = stats::BIC(mod_new),
          model_updated = NA
        )
        history <- tibble::add_row(history, history_row)
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(del), "` in backward step ",
          iter, ": the term could not be removed. See ertte_remove_term() ",
          "for details."
        ))
        next
      }
      p_val <- if (use_ic) NA_real_ else .ertte_anova_p(mod, mod_new)
      history_row <- tibble::tibble(
        iteration = iter,
        attempt = attm,
        step = "backward",
        criterion = criterion,
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
      if (use_ic) {
        candidate_ic <- as.numeric(ic_fn(mod_new))
        if (candidate_ic < best_metric) {
          update_ind <- attm
          best_metric <- candidate_ic
          best_mod <- mod_new
        }
      } else if (is.na(p_val)) {
        rlang::warn(paste0(
          "Skipping candidate term `", deparse(del), "` in backward step ",
          iter, ": comparison p-value is NA (often caused by a candidate ",
          "that's aliased/collinear with another term in the model)."
        ))
      } else if (p_val > best_metric) {
        update_ind <- attm
        best_metric <- p_val
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

# Refits `mod` with a new `formula`/`data`, dispatching on engine so
# `ertte_add_term()`/`ertte_remove_term()` (and therefore the SCM
# functions built on them) work across both `ertte_aft` and
# `ertte_coxph` models. `mod` itself is only used to determine which
# constructor to call (and, for AFT, which `dist` to refit with) --
# `stats::update()` doesn't work here, since the fitted object's `$call`
# refers to the constructor's own local argument bindings, not anything
# visible in the caller's frame.
.ertte_refit <- function(mod, formula, data) {
  UseMethod(".ertte_refit")
}

.ertte_refit.ertte_aft <- function(mod, formula, data) {
  ertte_aft(formula = formula, data = data, dist = mod$ertte$type)
}

.ertte_refit.ertte_coxph <- function(mod, formula, data) {
  ertte_coxph(formula = formula, data = data)
}

.ertte_refit.default <- function(mod, formula, data) {
  rlang::abort(paste0(
    "Don't know how to refit an object of class ", .fmt_bad_value(class(mod)),
    " -- ertte_add_term()/ertte_remove_term() currently only support ",
    "\"ertte_aft\" and \"ertte_coxph\" models."
  ))
}

#' Add or remove a covariate term from an exposure-response TTE model
#'
#' Add or remove a single covariate term from an existing ertte model,
#' returning a new fitted model object.
#'
#' @param mod An ertte model object.
#' @param term A one-sided formula naming the term to add/remove, e.g.
#' `~ sex`.
#' @param quiet Should the warning issued when the term can't be
#' added/removed (because it's already in the model / isn't in the model,
#' respectively) be suppressed? Defaults to `FALSE`.
#'
#' @details These functions are not typically called directly; they
#' underpin [ertte_scm_forward()] and [ertte_scm_backward()], used to
#' add or remove a single term from a time-to-event regression model.
#' `mod` is refit using the matching engine constructor based on its
#' class, so these functions work for both `ertte_aft` and `ertte_coxph`
#' models.
#' Regardless of whether the model is a parametric AFT model or a Cox
#' proportional hazard model, the `term` to be added or removed is 
#' defined by a single one-sided formula. Categorical covariates enter
#' the model as factor levels, whereas continuous covariates enter as 
#' linear terms by default. However, for a power-function parameterisation 
#' (`theta` such that `T ~ (x /ref)^theta` on the AFT time scale, or 
#' `h(t|x) ~ h0(t) * (x / ref)^theta` on the Cox hazard scale), the 
#' covariate can be wrapped in [ertte_power()], e.g. `~ ertte_power(age)`.
#'
#' @returns An ertte model object. If the term can't be added/removed
#' (see `quiet`), the original `mod` is returned unchanged.
#'
#' @name ertte_term
#' @examples
#' mod1 <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
#' mod2 <- ertte_add_term(mod1, ~ sex)
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
  .ertte_refit(mod, formula = fml, data = dat)
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
  # `survreg()` happens to accept a `terms` object directly as `formula`,
  # but `coxph()` doesn't (errors in `terms.formula()`/`ExtractVars`) --
  # convert to a plain formula so this works for both engines.
  .ertte_refit(mod, formula = stats::formula(trm_new), data = dat)
}
