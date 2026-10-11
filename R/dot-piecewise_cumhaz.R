#' Cumulative hazard of a piecewise exponential distribution (internal)
#'
#' @param x Numeric vector of non-negative times.
#' @param spec One element of the output of \code{.dropout_setup()}.
#' @return Numeric vector of cumulative hazards.
#' @keywords internal
#' @noRd
.piecewise_cumhaz <- function(x, spec) {
  if (length(spec$time) == 1L) return(spec$hazard * x)
  i <- findInterval(x, spec$time)
  spec$cum[i] + spec$hazard[i] * (x - spec$time[i])
}
