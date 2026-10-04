#' Response threshold of the latent normal variable (internal)
#'
#' Returns the value \code{c} such that \code{P(Z1 > z_tau, Z2 > c) = orr},
#' where \code{Z1} and \code{Z2} are standard normal with correlation
#' \code{theta}. Without a landmark (\code{z_tau = -Inf}) the threshold is
#' \code{qnorm(1 - orr)}.
#'
#' @param theta Copula correlation.
#' @param orr Objective response rate.
#' @param z_tau Latent PFS score of the landmark time.
#' @return Scalar threshold.
#' @keywords internal
#' @noRd
.resp_threshold <- function(theta, orr, z_tau) {
  if (z_tau == -Inf) {
    return(stats::qnorm(orr, lower.tail = FALSE))
  }
  f <- function(cc) .bvn_upper(z_tau, cc, theta) - orr
  stats::uniroot(f, interval = c(-12, 12), tol = 1e-12)$root
}
