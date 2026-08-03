test_that("er_predict.ertte_model() returns the same shape as ertte_predict()", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- er_predict.ertte_model(mod, ertte_data[1:5, ], time = 60)
  ref <- ertte_predict(mod, ertte_data[1:5, ], time = 60)
  # `expect_equal(ignore_attr = TRUE)` because row-subsetting a tibble
  # doesn't reliably preserve incidental per-column `label` attributes
  # (from `ertte_data`'s documentation metadata) -- irrelevant here.
  expect_equal(pred, ref, ignore_attr = TRUE)
})

test_that("er_predict.ertte_model() requires a time argument", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(er_predict.ertte_model(mod, ertte_data[1:5, ]), "time")
})

test_that("er_simulate.ertte_model() matches .ertte_simulate_draws()", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  sim <- er_simulate.ertte_model(mod, ertte_data, nsim = 5, seed = 99)
  expect_equal(sim, .ertte_simulate_draws(mod, ertte_data, nsim = 5, seed = 99), ignore_attr = TRUE)
})

test_that("er_summary.ertte_model() returns p_value/coefficients/glance", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  smm <- er_summary.ertte_model(mod)
  expect_named(smm, c("p_value", "coefficients", "glance"))
  expect_equal(nrow(smm$coefficients), 2L)
  expect_equal(nrow(smm$glance), 1L)
})

test_that("ertte registers er_predict/er_simulate/er_summary with erplots, if installed", {
  skip_if_not_installed("erplots")
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_true(!is.null(getS3method("er_predict", "ertte_model", envir = asNamespace("erplots"))))
})
