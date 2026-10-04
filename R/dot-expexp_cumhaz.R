#' Cumulative post-progression hazard of the exponential-exponential model
#' (internal)
#'
#' Returns a function that evaluates \code{G(t)} by cubic spline interpolation
#' of the tabulated values, with linear extension (slope \code{lam_o}) beyond
#' the grid. The spline keeps the integrands of the analytical functions
#' smooth, so that \code{integrate()} converges. The C++ generator
#' interpolates the same table linearly; the two differ by less than the error
#' of the tabulation itself.
#'
#' @param arm An \code{OncoArm} object with \code{os.model = "expexp"}.
#' @return A function of a numeric vector of times.
#' @keywords internal
#' @noRd
.expexp_cumhaz <- function(arm) {
  g <- arm$grid
  tmax <- g$t[length(g$t)]
  g_end <- g$G[length(g$G)]
  lam_o <- arm$lam_o
  sp <- stats::splinefun(g$t, g$G, method = "fmm")
  function(t) {
    out <- sp(pmin(t, tmax))
    over <- t > tmax
    out[over] <- g_end + lam_o * (t[over] - tmax)
    out
  }
}
