#' Calibrate the post-progression hazard to a target OS median (internal)
#'
#' Solves \code{S_OS(os.median) = 1/2} for the post-progression hazard of
#' non-responders in the illness-death model; responders have
#' \code{pps.hr.resp} times this hazard.
#'
#' @param arm A partially built \code{OncoArm} object.
#' @param os_median Target OS median.
#' @return Scalar hazard.
#' @keywords internal
#' @noRd
.calib_gam0 <- function(arm, os_median) {
  s_low <- 1 - arm$death.prop * (1 - exp(-arm$lam_p * os_median))
  if (s_low <= 0.5) {
    stop("'os.median' = ", os_median, " cannot be attained: even without ",
         "deaths after progression the OS survival at this time is below 0.5. ",
         "Reduce 'death.prop' or 'os.median'.", call. = FALSE)
  }
  f <- function(lg) .surv_os(os_median, arm, gam0 = exp(lg)) - 0.5
  exp(stats::uniroot(f, interval = c(log(1e-8), log(1e3)), tol = 1e-12)$root)
}
