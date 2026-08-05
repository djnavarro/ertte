test_that("ertte_aft_select_distribution() picks the true generating distribution", {
  # ertte_data's time/event were generated from a Weibull AFT model
  cmp <- ertte_aft_select_distribution(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_s3_class(cmp$comparison, "tbl_df")
  expect_setequal(cmp$comparison$dist, c("exponential", "weibull", "lognormal", "loglogistic"))
  expect_true(all(diff(cmp$comparison$aic) >= 0)) # sorted ascending by AIC
  expect_identical(cmp$comparison$dist[1], "weibull")
  expect_s3_class(cmp$model, "ertte_aft")
  expect_s3_class(cmp$model, "ertte_model")
  expect_identical(cmp$model$ertte$type, "weibull")
})

test_that("ertte_aft_select_distribution() validates candidates", {
  expect_error(
    ertte_aft_select_distribution(survival::Surv(time, event) ~ aucss, ertte_data, candidates = "bogus"),
    "`dist`"
  )
})

test_that("ertte_aft_select_distribution() rejects empty/missing candidates (issue #3)", {
  expect_error(
    ertte_aft_select_distribution(survival::Surv(time, event) ~ aucss, ertte_data, candidates = character(0)),
    "non-empty"
  )
  expect_error(
    ertte_aft_select_distribution(survival::Surv(time, event) ~ aucss, ertte_data, candidates = NULL),
    "non-empty"
  )
  expect_error(
    ertte_aft_select_distribution(survival::Surv(time, event) ~ aucss, ertte_data, candidates = c("weibull", NA)),
    "non-empty"
  )
})
