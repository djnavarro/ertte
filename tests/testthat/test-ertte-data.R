test_that("ertte_data has the expected columns and types", {
  expect_true(all(c(
    "id", "sex", "age", "weight", "dose", "treatment",
    "aucss", "cmaxss", "time", "event", "admin_censor"
  ) %in% names(ertte_data)))
  expect_true(all(ertte_data$time > 0))
  expect_true(all(ertte_data$event %in% c(0, 1)))
  expect_s3_class(ertte_data$sex, "factor")
  expect_s3_class(ertte_data$treatment, "factor")
  # fixed, known study-wide administrative cutoff -- same for every
  # subject, unlike `time` (see `.make_ertte_data()`)
  expect_true(all(ertte_data$admin_censor == 180))
  expect_true(all(ertte_data$time <= ertte_data$admin_censor + 1e-8))
})
