#' OS survival and density overall or by response (internal)
#'
#' For \code{os.model = "idm"} with response-independent post-progression
#' hazard and \code{response = "all"}, closed forms are used (Fleischer et al.
#' 2009, Theorem 5). Otherwise the contribution of patients who progressed is
#' integrated over the progression time \code{s}, split at the landmark.
#'
#' @param t Numeric vector of times.
#' @param arm An \code{OncoArm} object.
#' @param response One of \code{"all"}, \code{"responders"},
#'   \code{"nonresponders"}.
#' @param type \code{"survival"} or \code{"density"}.
#' @param gam0 Optional post-progression hazard of non-responders overriding
#'   the value stored in \code{arm} (used during calibration).
#' @return Numeric vector.
#' @keywords internal
#' @noRd
.surv_os <- function(t, arm, response = "all", type = "survival", gam0 = NULL) {
  lam_p <- arm$lam_p
  p <- arm$orr
  pi_d <- arm$death.prop
  tau <- arm$tau
  surv <- type == "survival"
  q_of_s <- function(s) .q_resp(.z_of_time(s, lam_p), arm$theta, arm$c_resp, arm$z_tau)

  if (arm$os.model == "idm") {
    g0 <- if (is.null(gam0)) arm$gam0 else gam0
    g1 <- arm$pps.hr.resp * g0
    if (response == "all" && arm$pps.hr.resp == 1) {
      if (abs(lam_p - g0) < 1e-10 * lam_p) {
        if (surv) return(exp(-lam_p * t) * (1 + (1 - pi_d) * lam_p * t))
        return(exp(-lam_p * t) * (pi_d * lam_p + (1 - pi_d) * lam_p * g0 * t))
      }
      if (surv) {
        return(exp(-lam_p * t) + (1 - pi_d) * lam_p *
                 (exp(-g0 * t) - exp(-lam_p * t)) / (lam_p - g0))
      }
      return(pi_d * lam_p * exp(-lam_p * t) + (1 - pi_d) * lam_p * g0 *
               (exp(-g0 * t) - exp(-lam_p * t)) / (lam_p - g0))
    }
    integrand <- function(s, tt) {
      q <- q_of_s(s)
      fp <- lam_p * exp(-lam_p * s)
      k0 <- if (surv) exp(-g0 * (tt - s)) else g0 * exp(-g0 * (tt - s))
      k1 <- if (surv) exp(-g1 * (tt - s)) else g1 * exp(-g1 * (tt - s))
      switch(response,
             all = fp * ((1 - q) * k0 + q * k1),
             responders = fp * q * k1 / p,
             nonresponders = fp * (1 - q) * k0 / (1 - p))
    }
  } else {
    if (response == "all") {
      if (surv) return(exp(-arm$lam_o * t))
      return(arm$lam_o * exp(-arm$lam_o * t))
    }
    cumhaz <- .expexp_cumhaz(arm)
    integrand <- function(s, tt) {
      q <- q_of_s(s)
      w <- if (response == "responders") q / p else (1 - q) / (1 - p)
      fp <- lam_p * exp(-lam_p * s)
      pis <- (arm$lam_o / lam_p) * exp(-arm$c_dec * s)
      k <- exp(-(cumhaz(tt) - cumhaz(s)))
      if (!surv) k <- k * .expexp_h12(tt, lam_p, arm$lam_o, arm$c_dec)
      fp * w * (1 - pis) * k
    }
  }

  # contribution of patients who were still progression-free (survival) or who
  # died before progression (density)
  first <- .surv_pfs(t, arm, response = response,
                     type = if (surv) "survival" else "density")
  if (!surv) {
    if (arm$os.model == "idm") {
      first <- pi_d * first
    } else {
      first <- (arm$lam_o / lam_p) * exp(-arm$c_dec * t) * first
    }
  }
  post <- vapply(t, function(tt) {
    if (tt <= 0) return(0)
    cuts <- if (tau > 0 && tau < tt) c(0, tau, tt) else c(0, tt)
    v <- 0
    for (i in seq_len(length(cuts) - 1L)) {
      v <- v + stats::integrate(integrand, lower = cuts[i], upper = cuts[i + 1L],
                                tt = tt, rel.tol = 1e-10, abs.tol = 1e-14,
                                subdivisions = 1000L)$value
    }
    v
  }, numeric(1))
  if (arm$os.model == "idm") post <- (1 - pi_d) * post
  first + post
}
