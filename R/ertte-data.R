

.make_ertte_data <- function(seed) {
  n <- 300L
  admin_censor <- 180 # days of administrative follow-up
  withr::with_seed(
    seed = seed,
    code = {
      ertte_data <- tibble::tibble(
        id = 1:n,
        sex = factor(sample(rep(c("Male", "Female"), c(n/2, n/2)))),
        age = sample(18:35, size = n, replace = TRUE)
      ) |>
        dplyr::mutate(
          weight = dplyr::if_else(
            condition = sex == "Male",
            true = (stats::runif(dplyr::n(), .05, .95)) |>
              stats::qlnorm(meanlog = 4.284, sdlog = 0.164) |>
              round(),
            false = (stats::runif(dplyr::n(), .05, .95)) |>
              stats::qlnorm(meanlog = 4.114, sdlog = 0.164) |>
              round()
          ),
          .by = sex
        ) |>
        dplyr::mutate(
          dose = sample(rep(c(0, 100, 200), c(n/3, n/3, n/3))),
          treatment = factor(dose == 0, levels = c(TRUE, FALSE), labels = c("Placebo", "Drug")),
          aucss = (stats::runif(n, .05, .95)) |>
            stats::qlnorm() |>
            (\(x) x * (dose + 10 * weight))() |>
            (\(x) dplyr::if_else(dose == 0, 0, x))() |>
            round(digits = 3),
          cmaxss = (exp(log(aucss/10) + stats::rnorm(n)/3) + stats::rnorm(n)) |>
            (\(x) dplyr::if_else(dose == 0, 0, x))() |>
            round(digits = 3),
        ) |>
        # simulated time-to-event outcome: a Weibull AFT model with an
        # exposure effect (higher aucss -> shorter time-to-event) and a
        # sex effect, plus administrative censoring at `admin_censor`
        # days and a little independent random (dropout) censoring.
        dplyr::mutate(
          mu = 5.2 - 0.0007 * aucss - 0.4 * as.numeric(sex == "Female"),
          scale = 0.7,
          time_true = exp(mu + scale * log(-log(1 - stats::runif(n)))),
          dropout_time = stats::rexp(n, rate = 1 / 400),
          time = pmax(pmin(time_true, dropout_time, admin_censor), 0.1) |> round(digits = 1),
          event = as.numeric(time_true <= pmin(dropout_time, admin_censor)),
        ) |>
        dplyr::select(-mu, -scale, -time_true, -dropout_time)
    }
  )
  attr(ertte_data$id, "label") <- "Subject"
  attr(ertte_data$sex, "label") <- "Sex"
  attr(ertte_data$age, "label") <- "Age"
  attr(ertte_data$weight, "label") <- "Weight"
  attr(ertte_data$dose, "label") <- "Dose"
  attr(ertte_data$treatment, "label") <- "Treatment"
  attr(ertte_data$aucss, "label") <- "AUCss"
  attr(ertte_data$cmaxss, "label") <- "Cmax,ss"
  attr(ertte_data$time, "label") <- "Time to event or censoring (days)"
  attr(ertte_data$event, "label") <- "Event indicator (1 = event, 0 = censored)"
  return(ertte_data)
}

#ertte_data <- .make_ertte_data(seed = 111L)
#usethis::use_data(ertte_data, overwrite = TRUE)

#' Sample simulated data for exposure-response time-to-event models
#'
#' @name ertte_data
#' @format A data frame with columns:
#' \describe{
#' \item{id}{Identifier}
#' \item{sex}{Sex}
#' \item{age}{Age}
#' \item{weight}{Weight}
#' \item{dose}{Nominal dose, units not specified}
#' \item{treatment}{Treatment}
#' \item{aucss}{AUCss}
#' \item{cmaxss}{Cmax,ss}
#' \item{time}{Time to event or right-censoring, in days}
#' \item{event}{Event indicator (1 = event observed, 0 = right-censored)}
#' }
#' @details
#'
#' This simulated dataset is entirely synthetic. `time`/`event` were
#' generated from a Weibull AFT model with an exposure (`aucss`) and sex
#' effect, administrative censoring at 180 days, and a little
#' independent random (dropout) censoring. You can find the data
#' generating code in the package source code.
#'
#' @examples
#' ertte_data
"ertte_data"
