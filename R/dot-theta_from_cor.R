#' Copula correlation that gives a target PFS and response correlation (internal)
#'
#' @param target Target \code{Corr(PFS, R)}.
#' @param orr Objective response rate.
#' @param z_tau Latent PFS score of the landmark time.
#' @return Scalar copula correlation in (-0.9999, 0.9999).
#' @keywords internal
#' @noRd
.theta_from_cor <- function(target, orr, z_tau) {
  f <- function(th) .cor_pfs_resp(th, orr, z_tau)$cor - target
  lim <- 0.9999
  f_lo <- f(-lim)
  f_hi <- f(lim)
  if (f_lo > 0 || f_hi < 0) {
    stop("'resp.cor' = ", target, " is outside the attainable range (",
         round(f_lo + target, 4), ", ", round(f_hi + target, 4),
         ") for this response rate and landmark. See CorBoundPFSResponse().",
         call. = FALSE)
  }
  stats::uniroot(f, interval = c(-lim, lim), f.lower = f_lo, f.upper = f_hi,
                 tol = 1e-10)$root
}
