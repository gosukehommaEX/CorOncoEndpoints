#' Distribution function of the accrual time (internal)
#'
#' @param x Numeric vector of times since the start of accrual.
#' @param acc Output of \code{.accrual_setup()}.
#' @return Probability that a patient has entered by time \code{x}.
#' @keywords internal
#' @noRd
.accrual_cdf <- function(x, acc) {
  if (length(acc$time) == 1L) return(as.numeric(x >= 0))
  stats::approx(acc$time, acc$cum, xout = x, yleft = 0, yright = 1)$y
}
