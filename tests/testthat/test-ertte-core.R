test_that("ertte_model() reproduces direct survreg() results", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  ref <- survival::survreg(survival::Surv(time, event) ~ aucss, ertte_data, dist = "weibull")
  expect_equal(unname(coef(mod)), unname(coef(ref)))
  expect_equal(as.numeric(stats::logLik(mod)), as.numeric(stats::logLik(ref)))
  expect_s3_class(mod, "ertte_model")
  expect_s3_class(mod, "survreg")
  expect_identical(mod$ertte$type, "weibull")
})

test_that("ertte_model() validates dist", {
  expect_error(
    ertte_model(survival::Surv(time, event) ~ aucss, ertte_data, dist = "not_a_dist"),
    "`dist`"
  )
})

test_that("ertte_predict() returns a tidy tibble with sane survival probabilities", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
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

test_that("ertte_predict() validates conf_level and time", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_predict(mod, ertte_data, time = 60, conf_level = 1.5), "`conf_level`")
  expect_error(ertte_predict(mod, ertte_data, time = -1), "`time`")
  expect_error(ertte_predict(mod, ertte_data, time = "a"), "`time`")
})

test_that("ertte_fun() reproduces the fitted model's own predictions", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_fun <- ertte_fun(mod)
  s1 <- unname(mod_fun(time = 60))
  s2 <- ertte_predict(mod, ertte_data, time = 60)$fit_survival
  expect_equal(s1, s2)
})

test_that("ertte_fun() checks param length", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_fun <- ertte_fun(mod)
  expect_error(mod_fun(param = c(1, 2, 3), time = 60), "length")
})
