#' Expected number of events by calendar time
#'
#' @description
#' Computes the expected number of observed PFS or OS events by given calendar
#' times, for groups described by \code{\link{OncoArm}} objects, with
#' piecewise-uniform accrual and exponential dropout.
#'
#' @details
#' For group \code{j} with \code{n_j} patients, event density \code{f_j},
#' dropout hazard \code{d_j} and accrual distribution function \code{A}, the
#' expected number of events by calendar time \code{T} is the integral of
#' \code{n_j f_j(x) exp(-d_j x) A(T - x)} over follow-up time \code{x} from 0
#' to \code{T}. The event densities are tabulated on a grid of 2001 points over
#' \code{[0, max(time)]} and the integral is computed by the midpoint rule.
#'
#' @param arms An \code{OncoArm} object or a list of them, one per group.
#' @param n Sample sizes, one per group.
#' @param a.time,a.rate Accrual as in \code{\link{rOncoEndpoints}}.
#' @param d.hazard Dropout hazard, a scalar or one value per group.
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
ExpectedEvents <- function(arms, n, a.time = NULL, a.rate = NULL, d.hazard = 0,
                           endpoint = c("os", "pfs"), time) {
  endpoint <- match.arg(endpoint)
  if (!is.numeric(time) || any(!is.finite(time)) || any(time <= 0)) {
    stop("'time' must be a numeric vector of positive calendar times.", call. = FALSE)
  }
  setup <- .design_setup(arms, n, a.time, a.rate, d.hazard, endpoint,
                         tmax = max(time))
  res <- t(vapply(time, function(tcal) {
    v <- .design_integrals(setup, tcal)
    c(v$events, v$by_group)
  }, numeric(1 + length(setup$arms))))
  out <- data.frame(time = time, res, check.names = FALSE)
  names(out) <- c("time", "events", names(setup$arms))
  out
}
