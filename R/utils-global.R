
# non-standard-evaluation variables used inside dplyr pipelines --
# declared here so `R CMD check` doesn't flag them as undefined globals.
utils::globalVariables(c(
  "aic", "aucss", "dat_id", "dose", "dropout_time", "iteration",
  "model_updated", "mu", "row_id", "sex", "sim_event", "sim_id",
  "sim_time", "time_true", "weight"
))

`%||%` <- function(x, y) {
  if (is.null(x)) return(y)
  x
}

.pick_seed <- function() {999 + sample.int(9000, size = 1L)}

# Renders a bad argument value for an error message without risking a
# second error from an unprintable/zero-length/weird-class input.
.fmt_bad_value <- function(x) {
  tryCatch(
    paste(deparse(x), collapse = " "),
    error = function(e) paste0("<", class(x)[1], ">")
  )
}
