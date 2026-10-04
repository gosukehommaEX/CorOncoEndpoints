#' Survival, density and hazard of an endpoint on a grid (internal)
#'
#' Used by the design functions. The grid starts at 0 and is uniform.
#'
#' @param arm An \code{OncoArm} object.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param tmax Upper end of the grid.
#' @param n_grid Number of grid points.
#' @return A list with elements \code{t}, \code{S}, \code{f} and \code{h}.
#' @keywords internal
#' @noRd
.endpoint_grid <- function(arm, endpoint, tmax, n_grid = 2001L) {
  tt <- seq(0, tmax, length.out = n_grid)
  if (endpoint == "pfs") {
    S <- exp(-arm$lam_p * tt)
    f <- arm$lam_p * S
  } else {
    os <- .surv_os_grid(tt, arm)
    S <- os$S
    f <- os$f
  }
  list(t = tt, S = S, f = f, h = f / S)
}
