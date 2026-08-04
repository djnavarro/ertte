test_that("ertte_coxph() reproduces direct coxph() results", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  ref <- survival::coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_equal(unname(coef(mod)), unname(coef(ref)))
  expect_equal(as.numeric(stats::logLik(mod)), as.numeric(stats::logLik(ref)))
  expect_s3_class(mod, "ertte_coxph")
  expect_s3_class(mod, "ertte_model")
  expect_s3_class(mod, "coxph")
  expect_identical(mod$ertte$type, "coxph")
})

test_that("ertte_coxph() stores the fitting data on the object", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_identical(mod$data, ertte_data)
})

test_that("ertte_predict()/ertte_fun() have no method yet for ertte_coxph", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_predict(mod, time = 60), "no applicable method")
  expect_error(ertte_fun(mod), "no applicable method")
})

test_that("ertte_add_term()/ertte_remove_term() work for ertte_coxph models", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  mod2 <- ertte_add_term(mod, ~sex)
  expect_s3_class(mod2, "ertte_coxph")
  expect_true("sex" %in% attr(stats::terms(mod2), "term.labels"))

  mod3 <- ertte_remove_term(mod2, ~sex)
  expect_s3_class(mod3, "ertte_coxph")
  expect_false("sex" %in% attr(stats::terms(mod3), "term.labels"))
})

test_that("ertte_scm_forward()/ertte_scm_backward() work for ertte_coxph models", {
  mod0 <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age", "weight"), seed = 4821)
  expect_s3_class(mod1, "ertte_coxph")
  # ertte_data's ground truth includes a sex effect
  expect_true("sex" %in% attr(stats::terms(mod1), "term.labels"))

  mod2 <- ertte_coxph(survival::Surv(time, event) ~ aucss + sex + dose, ertte_data)
  mod3 <- ertte_scm_backward(mod2, candidates = c("sex", "dose"), seed = 912)
  labs <- attr(stats::terms(mod3), "term.labels")
  expect_true("aucss" %in% labs)
  expect_false("dose" %in% labs) # not part of the ground truth model
})
