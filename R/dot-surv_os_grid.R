#' OS survival and density of all patients on a grid (internal)
#'
#' Fast version of \code{.surv_os()} for an increasing grid that starts at 0,
#' used by the design functions. For \code{os.model = "idm"} with a
#' response-dependent post-progression hazard, the contribution of patients
#' who progressed with response status \code{r} is
#' \code{J_r(t)}, the integral over \code{s < t} of
#' \code{f_P(s) P(R = r | PFS = s) exp(-gam_r (t - s))}. Because the kernel is
#' exponential, \code{J_r} is accumulated from one grid point to the next,
#' \code{J_r(t_i) = exp(-gam_r (t_i - t_(i-1))) J_r(t_(i-1))} plus the integral
#' over \code{(t_(i-1), t_i]}, which is computed by Gauss-Legendre quadrature.
#' The intervals are split at the landmark, and the first interval is
#' subdivided geometrically towards 0, where \code{P(R = r | PFS = s)} is not
#' smooth. Other cases use the closed forms of \code{.surv_os()}.
#'
#' @param tt Increasing numeric vector with \code{tt[1] = 0}.
#' @param arm An \code{OncoArm} object.
#' @param n_gl Number of Gauss-Legendre nodes per interval.
#' @return A list with elements \code{S} and \code{f}.
#' @keywords internal
#' @noRd
.surv_os_grid <- function(tt, arm, n_gl = 10L) {
  if (length(tt) < 2L || tt[1] != 0 || any(diff(tt) <= 0)) {
    stop("'tt' must be increasing and start at 0.", call. = FALSE)
  }
  if (arm$os.model != "idm" || arm$pps.hr.resp == 1) {
    return(list(S = .surv_os(tt, arm, type = "survival"),
                f = .surv_os(tt, arm, type = "density")))
  }
  lam_p <- arm$lam_p
  pi_d <- arm$death.prop
  gam <- c(arm$gam0, arm$gam1)
  t_end <- tt[length(tt)]
  extra <- tt[2] * 2 ^ -(1:20)
  if (arm$tau > 0 && arm$tau < t_end) extra <- c(extra, arm$tau)
  edges <- sort(unique(c(tt, extra)))
  lo <- edges[-length(edges)]
  hi <- edges[-1]
  half <- 0.5 * (hi - lo)
  gl <- .gauss_legendre(n_gl)
  s <- outer(0.5 * (lo + hi), rep(1, n_gl)) + outer(half, gl$x)
  wts <- outer(half, gl$w)
  fp <- lam_p * exp(-lam_p * s)
  q <- matrix(.q_resp(.z_of_time(s, lam_p), arm$theta, arm$c_resp, arm$z_tau),
              nrow = nrow(s))
  jj <- matrix(0, length(edges), 2L)
  for (r in 1:2) {
    a <- if (r == 1L) fp * (1 - q) else fp * q
    b <- rowSums(wts * a * exp(-gam[r] * (hi - s)))
    decay <- exp(-gam[r] * (hi - lo))
    for (i in seq_along(b)) jj[i + 1L, r] <- decay[i] * jj[i, r] + b[i]
  }
  idx <- match(tt, edges)
  s_pfs <- exp(-lam_p * tt)
  list(S = s_pfs + (1 - pi_d) * (jj[idx, 1] + jj[idx, 2]),
       f = pi_d * lam_p * s_pfs + (1 - pi_d) * (gam[1] * jj[idx, 1] + gam[2] * jj[idx, 2]))
}
