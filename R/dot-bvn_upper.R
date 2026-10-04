#' Upper orthant probability of a standard bivariate normal (internal)
#'
#' Computes \code{P(Z1 > a, Z2 > c)} for a standard bivariate normal vector
#' with correlation \code{theta}, by integrating the conditional probability of
#' \code{Z2} given \code{Z1} over \code{Z1 > a}. The integral is split at the
#' point where the conditional probability changes fastest.
#'
#' @param a,c Scalar thresholds (may be infinite).
#' @param theta Correlation in (-1, 1).
#' @return Scalar probability.
#' @keywords internal
#' @noRd
.bvn_upper <- function(a, c, theta) {
  if (a == Inf || c == Inf) return(0)
  if (a == -Inf) return(stats::pnorm(c, lower.tail = FALSE))
  if (c == -Inf) return(stats::pnorm(a, lower.tail = FALSE))
  if (theta == 0) {
    return(stats::pnorm(a, lower.tail = FALSE) * stats::pnorm(c, lower.tail = FALSE))
  }
  s <- sqrt(1 - theta ^ 2)
  f <- function(z) stats::dnorm(z) * stats::pnorm((theta * z - c) / s)
  zs <- c / theta
  if (is.finite(zs) && zs > a && abs(zs) < 8) {
    v <- stats::integrate(f, lower = a, upper = zs, rel.tol = 1e-10,
                          abs.tol = 1e-14, subdivisions = 1000L)$value +
      stats::integrate(f, lower = zs, upper = Inf, rel.tol = 1e-10,
                       abs.tol = 1e-14, subdivisions = 1000L)$value
  } else {
    v <- stats::integrate(f, lower = a, upper = Inf, rel.tol = 1e-10,
                          abs.tol = 1e-14, subdivisions = 1000L)$value
  }
  min(max(v, 0), 1)
}
