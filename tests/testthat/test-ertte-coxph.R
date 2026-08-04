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

test_that("ertte_fun() has no method yet for ertte_coxph", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_fun(mod), "no applicable method")
})

test_that("ertte_coxph() fits with model = TRUE (required by survfit())", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_false(is.null(mod$model))
})

test_that("ertte_predict() on ertte_coxph returns a tidy tibble with sane survival probabilities", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- ertte_predict(mod, ertte_data[1:5, ], time = c(30, 60, 90))
  expect_s3_class(pred, "tbl_df")
  expect_equal(nrow(pred), 15L)
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))
  expect_true(all(pred$ci_lower <= pred$fit_survival))
  expect_true(all(pred$ci_upper >= pred$fit_survival))
  # survival probability should decrease with time, for a fixed subject
  one_subject <- pred[pred$id == pred$id[1], ]
  expect_true(all(diff(one_subject$fit_survival) <= 0))
})

test_that("ertte_predict() on ertte_coxph handles a single-row newdata", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- ertte_predict(mod, ertte_data[1, , drop = FALSE], time = c(30, 60, 90))
  expect_equal(nrow(pred), 3L)
  expect_true(all(diff(pred$fit_survival) <= 0))
})

test_that("ertte_predict() on ertte_coxph extrapolates times beyond the observed range", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- ertte_predict(mod, ertte_data[1, , drop = FALSE], time = c(180, 10000))
  # extrapolated survival should be held constant (flat step-function
  # extension), not zero/NA
  expect_equal(pred$fit_survival[1], pred$fit_survival[2])
})

test_that("ertte_predict() on ertte_coxph validates conf_level and time", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_predict(mod, ertte_data, time = 60, conf_level = 1.5), "`conf_level`")
  expect_error(ertte_predict(mod, ertte_data, time = -1), "`time`")
  expect_error(ertte_predict(mod, ertte_data, time = "a"), "`time`")
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
