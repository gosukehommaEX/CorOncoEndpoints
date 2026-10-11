#' Expected number of events by calendar time
#'
#' @description
#' Computes the expected number of observed PFS or OS events by given calendar
#' times, for groups described by \code{\link{OncoArm}} objects, with
#' piecewise-uniform accrual, exponential or piecewise exponential dropout and,
#' for PFS, censoring of PFS.
#'
#' @details
#' For group \code{j} with \code{n_j} patients, event density \code{f_j},
#' cumulative dropout hazard \code{D_j} and accrual distribution function
#' \code{A}, the expected number of events by calendar time \code{T} is the
#' integral of \code{n_j f_j(x) exp(-D_j(x)) A(T - x)} over follow-up time
#' \code{x} from 0 to \code{T}; for PFS the integrand is also multiplied by
#' \code{exp(-c_j x)}, where \code{c_j} is \code{pfs.cens.hazard}. \code{A} is
#' piecewise linear, as specified by \code{a.time}, \code{a.rate} and
#' \code{a.prop} in \code{\link{rOncoEndpoints}}; the generated data have a
#' fixed number of patients in each accrual interval, which follows \code{A}
#' up to rounding. The event densities are tabulated on a grid of 2001 points
#' over \code{[0, max(time)]} and the integral is computed by the midpoint
#' rule.
#'
#' @param arms An \code{OncoArm} object or a list of them, one per group.
#' @param n Sample sizes, one per group.
#' @param a.time,a.rate,a.prop Accrual as in \code{\link{rOncoEndpoints}}.
#' @param d.hazard,d.median,d.time Dropout as in \code{\link{rOncoEndpoints}}.
#' @param pfs.cens.hazard Hazard of the censoring of PFS, a scalar or one
#'   value per group (used for \code{endpoint = "pfs"} only).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param time Numeric vector of calendar times.
#' @return A data frame with columns \code{time}, \code{events} (total) and one
#'   column per group.
#' @examples
#' ctl <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15, pps.hr.resp = 0.6)
#' trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.4,
#'                death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
#' ExpectedEvents(list(Control = ctl, Treatment = trt), n = c(350, 350),
#'                a.time = c(0, 24), endpoint = "os", time = c(12, 24, 36))
#' @export
ExpectedEvents <- function(arms, n, a.time = NULL, a.rate = NULL, a.prop = NULL,
                           d.hazard = 0, d.median = NULL, d.time = NULL,
                           pfs.cens.hazard = 0, endpoint = c("os", "pfs"), time) {
  endpoint <- match.arg(endpoint)
  if (!is.numeric(time) || any(!is.finite(time)) || any(time <= 0)) {
    stop("'time' must be a numeric vector of positive calendar times.", call. = FALSE)
  }
  setup <- .design_setup(arms, n, a.time, a.rate, a.prop, d.hazard, d.median,
                         d.time, pfs.cens.hazard, endpoint, tmax = max(time))
  res <- t(vapply(time, function(tcal) {
    v <- .design_integrals(setup, tcal)
    c(v$events, v$by_group)
  }, numeric(1 + length(setup$arms))))
  out <- data.frame(time = time, res, check.names = FALSE)
  names(out) <- c("time", "events", names(setup$arms))
  out
}
