test_that("ertte_scm_forward() adds the true generating term", {
  mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age", "weight"), seed = 4821)
  # ertte_data's ground truth includes a sex effect
  expect_true("sex" %in% attr(stats::terms(mod1), "term.labels"))
})

test_that("ertte_scm_backward() removes non-significant terms", {
  mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss + sex + dose, ertte_data)
  mod1 <- ertte_scm_backward(mod0, candidates = c("sex", "dose"), seed = 912)
  labs <- attr(stats::terms(mod1), "term.labels")
  expect_true("aucss" %in% labs)
  expect_false("dose" %in% labs) # not part of the ground truth model
})

test_that("ertte_scm_history() has the expected shape", {
  mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  h0 <- ertte_scm_history(mod0)
  expect_equal(nrow(h0), 1L)
  expect_equal(h0$iteration, 0L)

  mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"), seed = 55)
  h1 <- ertte_scm_history(mod1)
  expect_true(all(c(
    "iteration", "attempt", "step", "action", "term_tested", "model_tested",
    "model_converged", "term_p_value", "model_aic", "model_bic", "model_updated"
  ) %in% names(h1)))
})

test_that("SCM result doesn't depend on seed (only intermediate history row order might)", {
  mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_a <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age"), seed = 1)
  mod_b <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age"), seed = 2)
  expect_equal(
    sort(attr(stats::terms(mod_a), "term.labels")),
    sort(attr(stats::terms(mod_b), "term.labels"))
  )
})

test_that("ertte_add_term()/ertte_remove_term() validate `term`", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_add_term(mod, NULL), "`term`")
  expect_error(ertte_add_term(mod, time ~ sex), "two-sided")
  expect_error(ertte_add_term(mod, ~ sex + dose), "exactly one")
})

test_that("ertte_add_term()/ertte_remove_term() warn on no-ops", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_warning(ertte_add_term(mod, ~aucss), "already exists")
  expect_warning(ertte_remove_term(mod, ~sex), "does not exist")
  mod2 <- ertte_add_term(mod, ~sex, quiet = TRUE)
  expect_silent(ertte_remove_term(mod2, ~sex, quiet = TRUE))
})

test_that("ertte_scm_forward()/ertte_scm_backward() validate candidates", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_scm_forward(mod, candidates = "sex + dose"), "exactly one")
  expect_error(ertte_scm_backward(mod, candidates = character(0)), "non-empty")
})

test_that("ertte_scm_forward() skips (rather than misattributes) a candidate referencing a missing variable", {
  # Issue #8: a candidate whose variable isn't in the data used to make
  # ertte_add_term() silently no-op, leading to an anova() NA p-value
  # misreported as "aliased/collinear" rather than "variable not found".
  mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_warning(
    mod1 <- ertte_scm_forward(
      mod0,
      candidates = c("not_a_real_column", "sex"),
      seed = 77
    ),
    "not present in the fitting data"
  )
  # the valid candidate is still tried and selected
  expect_true("sex" %in% attr(stats::terms(mod1), "term.labels"))
  h1 <- ertte_scm_history(mod1)
  skipped_row <- h1[h1$term_tested %in% "~not_a_real_column", ]
  expect_true(nrow(skipped_row) >= 1L)
  expect_true(all(is.na(skipped_row$term_p_value)))
})

test_that("ertte_scm_forward()/ertte_scm_backward() skip a candidate whose refit errors, instead of crashing", {
  # Issue #9: a single-level factor candidate crashed the whole search
  # (contrasts error), rather than being skipped like other bad candidates.
  d <- ertte_data
  d$sex_single <- factor(rep("male", nrow(d)))
  mod0 <- ertte_aft(survival::Surv(time, event) ~ aucss, d)
  expect_warning(
    mod1 <- ertte_scm_forward(
      mod0,
      candidates = c("sex_single", "dose"),
      seed = 88
    ),
    "refitting the model with this term failed"
  )
  # the search completed and considered the other candidate
  h1 <- ertte_scm_history(mod1)
  expect_true("~dose" %in% h1$term_tested)
  failed_row <- h1[h1$term_tested %in% "~sex_single", ]
  expect_true(nrow(failed_row) >= 1L)
  expect_false(any(failed_row$model_converged))
})

test_that("ertte_scm_forward() skip-on-error behaviour also works for ertte_coxph()", {
  d <- ertte_data
  d$sex_single <- factor(rep("male", nrow(d)))
  mod0 <- ertte_coxph(survival::Surv(time, event) ~ aucss, d)
  expect_warning(
    mod1 <- ertte_scm_forward(
      mod0,
      candidates = c("sex_single", "dose"),
      seed = 90
    ),
    "refitting the model with this term failed"
  )
  h1 <- ertte_scm_history(mod1)
  expect_true("~dose" %in% h1$term_tested)
})
