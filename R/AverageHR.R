#' Average hazard ratio between two groups
#'
#' @description
#' Computes an average of the time-varying hazard ratio of the second group
#' versus the first, weighted by the expected number of events, for an
#' analysis at given calendar times.
#'
#' @details
#' The average hazard ratio at calendar time \code{T} is
#' \deqn{AHR(T) = exp( int log HR(x) w(x) dx / int w(x) dx ),}
#' where \code{HR(x)} is the ratio of the hazards of the second and first
#' groups at follow-up time \code{x} and \code{w(x)} is the expected number of
#' events at \code{x} in both groups among patients followed for at least
#' \code{x} by \code{T} (see \code{\link{ExpectedEvents}}). When the hazard
#' ratio is constant, \code{AHR(T)} equals it. With this average hazard ratio,
#' the Schoenfeld formula gives the number of events required by the log-rank
#' test (see \code{\link{RequiredEvents}}).
#'
#' @inheritParams ExpectedEvents
#' @param arms A list of two \code{OncoArm} objects, control first.
#' @param n Sample sizes of the two groups.
#' @return A data frame with columns \code{time}, \code{events} (expected,
#'   total) and \code{ahr}.
#' @examples
#' ctl <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15, pps.hr.resp = 0.6)
#' trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.4,
#'                death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
#' AverageHR(list(Control = ctl, Treatment = trt), n = c(350, 350),
#'           a.time = c(0, 24), endpoint = "os", time = c(24, 36))
#' @export
AverageHR <- function(arms, n, a.time = NULL, a.rate = NULL, d.hazard = 0,
                      endpoint = c("os", "pfs"), time) {
  endpoint <- match.arg(endpoint)
  if (!is.list(arms) || inherits(arms, "OncoArm") || length(arms) != 2L) {
    stop("'arms' must be a list of two OncoArm objects (control first).", call. = FALSE)
  }
  if (!is.numeric(time) || any(!is.finite(time)) || any(time <= 0)) {
    stop("'time' must be a numeric vector of positive calendar times.", call. = FALSE)
  }
  setup <- .design_setup(arms, n, a.time, a.rate, d.hazard, endpoint,
                         tmax = max(time))
  res <- t(vapply(time, function(tcal) {
    v <- .design_integrals(setup, tcal)
    c(v$events, v$ahr)
  }, numeric(2)))
  data.frame(time = time, events = res[, 1], ahr = res[, 2])
}
