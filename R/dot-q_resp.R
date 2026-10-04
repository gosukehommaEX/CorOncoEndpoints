#' Conditional response probability given the latent PFS score (internal)
#'
#' Returns \code{P(R = 1 | Z1 = z)}, which is zero below the landmark score
#' \code{z_tau} and \code{pnorm((theta * z - c_resp) / sqrt(1 - theta^2))}
#' above it.
#'
#' @param z Numeric vector of latent PFS scores.
#' @param theta Copula correlation.
#' @param c_resp Response threshold.
#' @param z_tau Latent PFS score of the landmark time.
#' @return Numeric vector of probabilities.
#' @keywords internal
#' @noRd
.q_resp <- function(z, theta, c_resp, z_tau) {
  s <- sqrt(1 - theta ^ 2)
  ifelse(z > z_tau, stats::pnorm((theta * z - c_resp) / s), 0)
}
