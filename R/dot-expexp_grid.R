#' Tabulated cumulative post-progression hazard of the exponential-exponential
#' model (internal)
#'
#' Tabulates \code{G(t)}, the integral of \code{.expexp_h12()} from 0 to
#' \code{t}, on a uniform grid over \code{[0, 40 / lam_o]} by the trapezoidal
#' rule. Beyond the grid, \code{G} increases linearly with slope \code{lam_o},
#' the limit of the hazard. The grid also stores the expected post-progression
#' survival \code{m(s)} for a progression at time \code{s}.
#'
#' @param lam_p,lam_o,c_dec Model parameters.
#' @param n_grid Number of grid points.
#' @return A list with elements \code{h} (grid step), \code{t}, \code{G} and
#'   \code{m}.
#' @keywords internal
#' @noRd
.expexp_grid <- function(lam_p, lam_o, c_dec, n_grid = 20001L) {
  tmax <- 40 / lam_o
  tt <- seq(0, tmax, length.out = n_grid)
  hz <- .expexp_h12(tt, lam_p, lam_o, c_dec)
  step <- tt[2] - tt[1]
  G <- c(0, cumsum(0.5 * (hz[-1] + hz[-n_grid]) * step))
  eG <- exp(-G)
  tail_int <- rev(cumsum(rev(c(0.5 * (eG[-1] + eG[-n_grid]) * step, 0))))
  tail_int <- tail_int + eG[n_grid] / lam_o
  m <- exp(G) * tail_int
  list(h = step, t = tt, G = G, m = m)
}
