#' Copula correlation that gives a target median PFS among responders (internal)
#'
#' @param m1 Target median PFS among responders.
#' @param orr Objective response rate.
#' @param lam_p PFS hazard.
#' @param z_tau Latent PFS score of the landmark time.
#' @return Scalar copula correlation in (-0.9999, 0.9999).
#' @keywords internal
#' @noRd
.theta_from_resp_median <- function(m1, orr, lam_p, z_tau) {
  z_m <- max(.z_of_time(m1, lam_p), z_tau)
  f <- function(th) {
    cc <- .resp_threshold(th, orr, z_tau)
    .bvn_upper(z_m, cc, th) / orr - 0.5
  }
  lim <- 0.9999
  f_lo <- f(-lim)
  f_hi <- f(lim)
  if (f_lo > 0 || f_hi < 0) {
    stop("'resp.pfs.median' = ", m1, " cannot be attained for this response ",
         "rate, PFS median and landmark.", call. = FALSE)
  }
  stats::uniroot(f, interval = c(-lim, lim), f.lower = f_lo, f.upper = f_hi,
                 tol = 1e-10)$root
}
