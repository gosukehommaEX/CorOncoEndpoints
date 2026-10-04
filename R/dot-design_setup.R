#' Common setup of the design functions (internal)
#'
#' Checks the inputs of \code{ExpectedEvents()}, \code{AverageHR()} and
#' \code{RequiredEvents()} and tabulates the event density and hazard of each
#' group on a uniform grid over \code{[0, tmax]}.
#'
#' @param arms,n,a.time,a.rate,d.hazard,endpoint See \code{ExpectedEvents()}.
#' @param tmax Upper end of the grid.
#' @param n_grid Number of grid points.
#' @return A list with elements \code{arms}, \code{n}, \code{d}, \code{acc}
#'   and \code{grids}.
#' @keywords internal
#' @noRd
.design_setup <- function(arms, n, a.time, a.rate, d.hazard, endpoint, tmax,
                          n_grid = 2001L) {
  arms <- .check_arms(arms)
  k <- length(arms)
  if (!is.numeric(n) || length(n) != k || any(n <= 0)) {
    stop("'n' must be a vector of positive sample sizes, one per group.", call. = FALSE)
  }
  if (!is.numeric(d.hazard) || !(length(d.hazard) %in% c(1L, k)) ||
      any(!is.finite(d.hazard)) || any(d.hazard < 0)) {
    stop("'d.hazard' must be a non-negative number or one per group.", call. = FALSE)
  }
  acc <- .accrual_setup(a.time, a.rate)
  grids <- lapply(arms, .endpoint_grid, endpoint = endpoint, tmax = tmax,
                  n_grid = n_grid)
  list(arms = arms, n = n, d = rep_len(d.hazard, k), acc = acc, grids = grids)
}
