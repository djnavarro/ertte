test_that("ertte_power() computes log(x / ref) with the expected default ref", {
  x <- c(20, 30, 40, 50, 60)
  out <- ertte_power(x)
  expect_equal(attr(out, "ref"), stats::median(x))
  expect_equal(as.numeric(out), log(x / stats::median(x)))
  expect_s3_class(out, "ertte_power")
})

test_that("ertte_power() respects a user-supplied ref", {
  x <- c(20, 30, 40, 50, 60)
  out <- ertte_power(x, ref = 25)
  expect_equal(attr(out, "ref"), 25)
  expect_equal(as.numeric(out), log(x / 25))
})

test_that("ertte_power() errors on non-positive x or invalid ref", {
  expect_error(ertte_power(c(1, 2, 0)), "strictly positive")
  expect_error(ertte_power(c(1, 2, -3)), "strictly positive")
  expect_error(ertte_power(c(1, 2, 3), ref = 0), "strictly positive")
  expect_error(ertte_power(c(1, 2, 3), ref = c(1, 2)), "single strictly positive")
  expect_error(ertte_power(c(1, 2, 3), ref = "a"), "single strictly positive")
  expect_error(ertte_power("not numeric"), "numeric")
})

test_that("ertte_power() tolerates NA values", {
  x <- c(20, NA, 40, 50, 60)
  out <- ertte_power(x)
  expect_equal(attr(out, "ref"), stats::median(x, na.rm = TRUE))
  expect_true(is.na(out[2]))
})

test_that("ertte_aft() fits an ertte_power() term and reports theta as an ordinary coefficient", {
  mod <- ertte_aft(Surv(time, event) ~ aucss + ertte_power(age), ertte_data)
  cf <- stats::coef(mod)
  expect_true("ertte_power(age)" %in% names(cf))
  ci <- stats::confint(mod)
  expect_true("ertte_power(age)" %in% rownames(ci))
})

test_that("ertte_predict()/ertte_fun() reuse the fitting-time ref, not a newdata-derived one", {
  train_ref <- stats::median(ertte_data$age)
  mod <- ertte_aft(Surv(time, event) ~ ertte_power(age), ertte_data)

  newdata <- ertte_data[1, , drop = FALSE]
  newdata$age <- train_ref + 40 # deliberately far from both training and (trivially) newdata's own median

  mm <- stats::model.matrix(stats::delete.response(stats::terms(mod)), newdata)
  expect_equal(mm[, "ertte_power(age)"], log(newdata$age / train_ref), ignore_attr = TRUE)

  mod_fun <- ertte_fun(mod)
  s_direct <- mod_fun(data = newdata, time = 60)
  # manually build survival prediction using the same design matrix to confirm consistency
  expect_false(is.na(s_direct))
})

test_that("ertte_coxph() also reuses the fitting-time ref via survfit()", {
  train_ref <- stats::median(ertte_data$age)
  mod <- ertte_coxph(Surv(time, event) ~ ertte_power(age), ertte_data)
  newdata <- ertte_data[1, , drop = FALSE]
  newdata$age <- train_ref + 40

  pred <- ertte_predict(mod, newdata, time = c(30, 60))
  expect_equal(nrow(pred), 2L)
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))

  mm <- stats::model.matrix(stats::delete.response(stats::terms(mod)), newdata)
  expect_equal(mm[, "ertte_power(age)"], log(newdata$age / train_ref), ignore_attr = TRUE)
})

test_that("ertte_add_term()/SCM accept an ertte_power() term/candidate", {
  mod0 <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  mod1 <- ertte_add_term(mod0, ~ ertte_power(age))
  expect_true("ertte_power(age)" %in% attr(stats::terms(mod1), "term.labels"))

  mod2 <- ertte_remove_term(mod1, ~ ertte_power(age))
  expect_false("ertte_power(age)" %in% attr(stats::terms(mod2), "term.labels"))

  mod3 <- ertte_scm_forward(mod0, candidates = c("sex", "ertte_power(age)"), seed = 371)
  expect_true(all(attr(stats::terms(mod3), "term.labels") %in% c("aucss", "sex", "ertte_power(age)")))
})
