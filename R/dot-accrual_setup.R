#' Check and tabulate a piecewise-uniform accrual distribution (internal)
#'
#' @param a.time Accrual breakpoints starting at 0, or \code{NULL} when all
#'   patients enter at time 0.
#' @param a.rate Relative accrual intensities on the intervals between the
#'   breakpoints (default: equal intensities).
#' @return A list with elements \code{time} (breakpoints) and \code{cum}
#'   (cumulative probabilities at the breakpoints).
#' @keywords internal
#' @noRd
.accrual_setup <- function(a.time, a.rate) {
  if (is.null(a.time)) {
    if (!is.null(a.rate)) stop("'a.rate' requires 'a.time'.", call. = FALSE)
    return(list(time = 0, cum = 0))
  }
  if (!is.numeric(a.time) || length(a.time) < 2L || a.time[1] != 0 ||
      any(diff(a.time) <= 0) || any(!is.finite(a.time))) {
    stop("'a.time' must be a strictly increasing finite numeric vector ",
         "starting at 0 with at least two elements.", call. = FALSE)
  }
  k <- length(a.time) - 1L
  if (is.null(a.rate)) a.rate <- rep(1, k)
  if (!is.numeric(a.rate) || length(a.rate) != k || any(a.rate < 0) ||
      sum(a.rate) <= 0) {
    stop("'a.rate' must have length(a.time) - 1 non-negative elements with a ",
         "positive sum.", call. = FALSE)
  }
  w <- a.rate * diff(a.time)
  list(time = a.time, cum = c(0, cumsum(w) / sum(w)))
}
