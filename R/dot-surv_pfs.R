#' PFS survival, density and hazard overall or by response (internal)
#'
#' @param t Numeric vector of times.
#' @param arm An \code{OncoArm} object.
#' @param response One of \code{"all"}, \code{"responders"},
#'   \code{"nonresponders"}.
#' @param type One of \code{"survival"}, \code{"density"}, \code{"hazard"}.
#' @return Numeric vector.
#' @keywords internal
#' @noRd
.surv_pfs <- function(t, arm, response = "all", type = "survival") {
  lam_p <- arm$lam_p
  p <- arm$orr
  s_all <- exp(-lam_p * t)
  f_all <- lam_p * s_all
  if (response == "all") {
    s <- s_all
    f <- f_all
  } else {
    s1 <- vapply(t, function(tt) {
      .bvn_upper(max(.z_of_time(tt, lam_p), arm$z_tau), arm$c_resp, arm$theta) / p
    }, numeric(1))
    q <- .q_resp(.z_of_time(t, lam_p), arm$theta, arm$c_resp, arm$z_tau)
    if (response == "responders") {
      s <- s1
      f <- f_all * q / p
    } else {
      s <- (s_all - p * s1) / (1 - p)
      f <- f_all * (1 - q) / (1 - p)
    }
  }
  switch(type, survival = s, density = f, hazard = f / s)
}
