test_that("simulate.ertte_model() has the expected shape", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  sim <- simulate(mod, nsim = 5, seed = 321)
  expect_s3_class(sim, "tbl_df")
  expect_equal(nrow(sim), nrow(ertte_data) * 5L)
  expect_true(all(c("dat_id", "sim_id", "sim_time", "sim_event", "coef_(Intercept)", "coef_aucss") %in% names(sim)))
  expect_true(all(sim$sim_event %in% c(0, 1)))
  expect_true(all(sim$sim_time <= ertte_data$time[sim$dat_id] + 1e-8))
})

test_that("simulate.ertte_model() is reproducible given a seed", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  sim_a <- simulate(mod, nsim = 5, seed = 42)
  sim_b <- simulate(mod, nsim = 5, seed = 42)
  expect_equal(sim_a, sim_b)
})

test_that("simulate.ertte_model() reports an auto-picked seed", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_message(simulate(mod, nsim = 5), "Using seed")
})

test_that(".ertte_simulate_draws() requires response columns in newdata", {
  mod <- ertte_model(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(
    .ertte_simulate_draws(mod, newdata = ertte_data[c("aucss")], nsim = 5, seed = 1),
    "must contain"
  )
})
