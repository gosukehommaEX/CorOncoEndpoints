#' Standard normal score of the PFS distribution function (internal)
#'
#' Returns \code{qnorm(1 - exp(-lam_p * t))}, the value of the latent normal
#' variable that corresponds to a progression-free survival time \code{t}.
#'
#' @param t Numeric vector of times (non-negative).
#' @param lam_p PFS hazard.
#' @return Numeric vector; \code{-Inf} at \code{t = 0}.
#' @keywords internal
#' @noRd
.z_of_time <- function(t, lam_p) {
  stats::qnorm(-expm1(-lam_p * t))
}
