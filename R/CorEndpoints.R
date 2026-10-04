#' Pearson correlations among PFS, OS and response
#'
#' @description
#' Computes the Pearson correlations \code{Corr(PFS, R)}, \code{Corr(OS, R)}
#' and \code{Corr(PFS, OS)} implied by an \code{OncoArm} object, together with
#' the partial correlation of OS and response given PFS.
#'
#' @details
#' Write \code{OS = PFS + Delta * D}, where \code{Delta} indicates progression
#' as the first event and \code{D} is the post-progression survival. For
#' \code{os.model = "idm"}, \code{Delta} and \code{D} are independent of PFS
#' given response, and with \code{m_r = (1 - death.prop) / gam_r} and
#' \code{delta = m_1 - m_0},
#' \deqn{Cov(OS, R) = Cov(PFS, R) + p (1 - p) delta,}
#' \deqn{Cov(PFS, OS) = Var(PFS) + delta Cov(PFS, R),}
#' \deqn{Var(OS) = Var(PFS) + Var(Delta D) + 2 delta Cov(PFS, R),}
#' where \code{p} is the response rate. The partial covariance of OS and
#' response given PFS is \code{delta p (1 - p) (1 - Corr(PFS, R)^2)}, so the
#' partial correlation has the sign of \code{delta}; it is zero when
#' \code{pps.hr.resp = 1}, in which case
#' \code{Corr(OS, R) = Corr(PFS, R) Corr(PFS, OS)}.
#' \code{Cov(PFS, R)} is a one-dimensional integral over the latent normal
#' score of PFS.
#'
#' For \code{os.model = "expexp"}, the variances are those of the exponential
#' distributions, and the covariances are computed from the expected
#' post-progression survival given the progression time, which is tabulated.
#'
#' @param arm An object of class \code{OncoArm}.
#' @return A named numeric vector with elements \code{cor.pfs.resp},
#'   \code{cor.os.resp}, \code{cor.pfs.os} and \code{pcor.os.resp} (partial
#'   correlation of OS and response given PFS).
#' @references
#' Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A statistical
#' model for the dependence between progression-free survival and overall
#' survival. \emph{Statistics in Medicine}, 28, 2669--2686.
#' \doi{10.1002/sim.3637}
#' @examples
#' arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15, pps.hr.resp = 0.6)
#' CorEndpoints(arm)
#' @export
CorEndpoints <- function(arm) {
  if (!inherits(arm, "OncoArm")) stop("'arm' must be an OncoArm object.", call. = FALSE)
  p <- arm$orr
  lam_p <- arm$lam_p
  rho_pr <- .cor_pfs_resp(arm$theta, p, arm$z_tau)$cor
  cov_pr <- rho_pr * sqrt(p * (1 - p)) / lam_p
  var_p <- 1 / lam_p ^ 2
  if (arm$os.model == "idm") {
    pd <- arm$death.prop
    m0 <- (1 - pd) / arm$gam0
    m1 <- (1 - pd) / arm$gam1
    s0 <- 2 * (1 - pd) / arm$gam0 ^ 2
    s1 <- 2 * (1 - pd) / arm$gam1 ^ 2
    dm <- m1 - m0
    var_dd <- p * s1 + (1 - p) * s0 - (p * m1 + (1 - p) * m0) ^ 2
    cov_or <- cov_pr + p * (1 - p) * dm
    cov_po <- var_p + dm * cov_pr
    var_o <- var_p + var_dd + 2 * dm * cov_pr
  } else {
    g <- arm$grid
    pis <- (arm$lam_o / lam_p) * exp(-arm$c_dec * g$t)
    # cubic spline keeps the integrands smooth (see .expexp_cumhaz())
    gfun <- stats::splinefun(g$t, (1 - pis) * g$m, method = "fmm")
    top <- g$t[length(g$t)]
    cuts <- if (arm$tau > 0) c(0, arm$tau, top) else c(0, top)
    int <- function(fn) {
      v <- 0
      for (i in seq_len(length(cuts) - 1L)) {
        v <- v + stats::integrate(fn, lower = cuts[i], upper = cuts[i + 1L],
                                  rel.tol = 1e-10, abs.tol = 1e-12,
                                  subdivisions = 2000L)$value
      }
      v
    }
    fp <- function(s) lam_p * exp(-lam_p * s)
    qs <- function(s) .q_resp(.z_of_time(s, lam_p), arm$theta, arm$c_resp, arm$z_tau)
    e_pdd <- int(function(s) s * fp(s) * gfun(s))
    e_dd <- int(function(s) fp(s) * gfun(s))
    cov_po <- var_p + e_pdd - e_dd / lam_p
    cov_or <- cov_pr + int(function(s) gfun(s) * (qs(s) - p) * fp(s))
    var_o <- 1 / arm$lam_o ^ 2
  }
  var_r <- p * (1 - p)
  r_or <- cov_or / sqrt(var_o * var_r)
  r_po <- cov_po / sqrt(var_p * var_o)
  pcor <- (r_or - rho_pr * r_po) / sqrt((1 - rho_pr ^ 2) * (1 - r_po ^ 2))
  c(cor.pfs.resp = rho_pr, cor.os.resp = r_or, cor.pfs.os = r_po,
    pcor.os.resp = pcor)
}
