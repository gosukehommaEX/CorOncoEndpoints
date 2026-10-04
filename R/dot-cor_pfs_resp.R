#' Correlation between PFS and response for a given copula correlation (internal)
#'
#' Computes \code{Corr(PFS, R)} for exponential PFS. The value does not depend
#' on the PFS hazard. With \code{lam_p * PFS = -log(Phi(-Z1))},
#' \code{Cov(lam_p * PFS, R)} is the integral of
#' \code{(-log(Phi(-z)) - 1) * P(R = 1 | Z1 = z) * dnorm(z)} over \code{z}.
#'
#' @param theta Copula correlation.
#' @param orr Objective response rate.
#' @param z_tau Latent PFS score of the landmark time.
#' @return A list with elements \code{cor} and \code{c_resp}.
#' @keywords internal
#' @noRd
.cor_pfs_resp <- function(theta, orr, z_tau) {
  c_resp <- .resp_threshold(theta, orr, z_tau)
  s <- sqrt(1 - theta ^ 2)
  g <- function(z) {
    (-stats::pnorm(-z, log.p = TRUE) - 1) *
      stats::pnorm((theta * z - c_resp) / s) * stats::dnorm(z)
  }
  lo <- z_tau
  cuts <- lo
  if (theta != 0) {
    zs <- c_resp / theta
    # split only where the normal density is not negligible: a split point far
    # in the tail makes integrate() miss the mass near zero
    if (zs > lo && abs(zs) < 8) cuts <- c(cuts, zs)
  }
  cuts <- c(cuts, Inf)
  v <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    v <- v + stats::integrate(g, lower = cuts[i], upper = cuts[i + 1L],
                              rel.tol = 1e-10, abs.tol = 1e-14,
                              subdivisions = 1000L)$value
  }
  list(cor = v / sqrt(orr * (1 - orr)), c_resp = c_resp)
}
