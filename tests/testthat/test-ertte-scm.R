test_that("ertte_scm_forward() adds the true generating term", {
  mod0 <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age", "weight"), seed = 4821)
  # ertte_data's ground truth includes a sex effect
  expect_true("sex" %in% attr(stats::terms(mod1), "term.labels"))
})

test_that("ertte_scm_backward() removes non-significant terms", {
  mod0 <- ertte_model(survival::Surv(time, event) ~ aucss + sex + dose, ertte_data)
  mod1 <- ertte_scm_backward(mod0, candidates = c("sex", "dose"), seed = 912)
  labs <- attr(stats::terms(mod1), "term.labels")
  expect_true("aucss" %in% labs)
  expect_false("dose" %in% labs) # not part of the ground truth model
})

test_that("ertte_scm_history() has the expected shape", {
  mod0 <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
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
  mod0 <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_a <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age"), seed = 1)
  mod_b <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age"), seed = 2)
  expect_equal(
    sort(attr(stats::terms(mod_a), "term.labels")),
    sort(attr(stats::terms(mod_b), "term.labels"))
  )
})

test_that("ertte_add_term()/ertte_remove_term() validate `term`", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_add_term(mod, NULL), "`term`")
  expect_error(ertte_add_term(mod, time ~ sex), "two-sided")
  expect_error(ertte_add_term(mod, ~ sex + dose), "exactly one")
})

test_that("ertte_add_term()/ertte_remove_term() warn on no-ops", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_warning(ertte_add_term(mod, ~aucss), "already exists")
  expect_warning(ertte_remove_term(mod, ~sex), "does not exist")
  mod2 <- ertte_add_term(mod, ~sex, quiet = TRUE)
  expect_silent(ertte_remove_term(mod2, ~sex, quiet = TRUE))
})

test_that("ertte_scm_forward()/ertte_scm_backward() validate candidates", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_scm_forward(mod, candidates = "sex + dose"), "exactly one")
  expect_error(ertte_scm_backward(mod, candidates = character(0)), "non-empty")
})
