#' Post-progression hazard of the exponential-exponential model (internal)
#'
#' With PFS hazard \code{lam_p}, OS hazard \code{lam_o} and pre-progression
#' death hazard \code{lam_o * exp(-c_dec * t)}, OS is exponential if and only
#' if the post-progression hazard (Markov, time since randomization) is
#' \code{lam_o * (1 - exp(-k t)) / (1 - exp(-b t))} with
#' \code{k = lam_p + c_dec - lam_o} and \code{b = lam_p - lam_o}.
#'
#' @param t Numeric vector of times.
#' @param lam_p,lam_o,c_dec Model parameters.
#' @return Numeric vector of hazards; the limit \code{lam_o * k / b} at zero.
#' @keywords internal
#' @noRd
.expexp_h12 <- function(t, lam_p, lam_o, c_dec) {
  k <- lam_p + c_dec - lam_o
  b <- lam_p - lam_o
  out <- lam_o * (-expm1(-k * t)) / (-expm1(-b * t))
  out[t < 1e-12] <- lam_o * k / b
  out
}
