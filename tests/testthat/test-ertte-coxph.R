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
