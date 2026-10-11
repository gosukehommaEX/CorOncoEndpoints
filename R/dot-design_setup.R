#' Common setup of the design functions (internal)
#'
#' Checks the inputs of \code{ExpectedEvents()}, \code{AverageHR()} and
#' \code{RequiredEvents()} and tabulates the event density and hazard of each
#' group on a uniform grid over \code{[0, tmax]}.
#'
#' @param arms,n,a.time,a.rate,a.prop,d.hazard,d.median,d.time,pfs.cens.hazard,endpoint
#'   See \code{ExpectedEvents()}.
#' @param tmax Upper end of the grid.
#' @param n_grid Number of grid points.
#' @return A list with elements \code{arms}, \code{n}, \code{drop} (dropout
#'   distributions), \code{pc} (hazards of the censoring of the endpoint other
#'   than dropout: those of PFS for \code{endpoint = "pfs"}, 0 for OS),
#'   \code{acc} and \code{grids}.
#' @keywords internal
#' @noRd
.design_setup <- function(arms, n, a.time, a.rate, a.prop = NULL, d.hazard,
                          d.median = NULL, d.time = NULL, pfs.cens.hazard = 0,
                          endpoint, tmax, n_grid = 2001L) {
  arms <- .check_arms(arms)
  k <- length(arms)
  if (!is.numeric(n) || length(n) != k || any(n <= 0)) {
    stop("'n' must be a vector of positive sample sizes, one per group.", call. = FALSE)
  }
  for (a in arms) .check_endpoint_available(a, endpoint)
  acc <- .accrual_setup(a.time, a.rate, a.prop, sum(n))
  drop <- .dropout_setup(d.hazard, d.median, d.time, k)
  pc <- .group_rates(pfs.cens.hazard, k, "pfs.cens.hazard")
  grids <- lapply(arms, .endpoint_grid, endpoint = endpoint, tmax = tmax,
                  n_grid = n_grid)
  list(arms = arms, n = n, drop = drop,
       pc = if (endpoint == "pfs") pc else rep(0, k), acc = acc, grids = grids)
}
