#' Objects exported from other packages
#'
#' These objects are imported from other packages. Follow the link below to
#' see its documentation.
#'
#' \describe{
#'   \item{survival}{\code{\link[survival]{Surv}}}
#' }
#'
#' @name reexports
#' @keywords internal
#' @return `Surv()` returns an object of class `"Surv"`, a matrix of
#'   time/event (or time/time2/event, for interval- or counting-process-style
#'   data) values used as the response in [ertte_aft()]/[ertte_coxph()]
#'   formulas. See `survival::Surv()` for full details on its return value.
#' @importFrom survival Surv
#' @export
survival::Surv
NULL
