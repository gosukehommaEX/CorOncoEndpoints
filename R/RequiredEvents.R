#' Number of events required by the log-rank test
#'
#' @description
#' Computes the number of PFS or OS events required for a one-sided log-rank
#' test to have a given power, using the Schoenfeld formula with the average
#' hazard ratio of \code{\link{AverageHR}}, and the calendar time at which this
#' number of events is expected.
#'
#' @details
#' With allocation proportions \code{r1} and \code{r2}, the required number of
#' events is
#' \deqn{D = (z_{1 - alpha} + z_{power})^2 / (r1 r2 log(AHR)^2).}
#' Because the average hazard ratio depends on the analysis time when the hazard
#' ratio varies over time (as it does for OS under \code{os.model = "idm"}),
#' the number of events and the analysis time are found by iteration: the
#' analysis time is the calendar time at which \code{D} events are expected,
#' and the average hazard ratio is recomputed at that time until the rounded
#' number of events no longer changes. For exponential PFS in both groups the
#' hazard ratio is constant and the result is the usual Schoenfeld number.
#'
#' @inheritParams AverageHR
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @param tmax Upper limit of calendar time searched (default: end of accrual
#'   plus 10 times the largest median of the endpoint).
#' @return A list with elements \code{events} (rounded up), \code{events.exact},
#'   \code{time} (calendar time at which \code{events} are expected),
#'   \code{ahr}, \code{max.events} (expected events by \code{tmax}) and
#'   \code{iterations}.
#' @examples
#' ctl <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15, pps.hr.resp = 0.6)
#' trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.4,
#'                death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
#' RequiredEvents(list(Control = ctl, Treatment = trt), n = c(350, 350),
#'                a.time = c(0, 24), endpoint = "os")
#' @export
RequiredEvents <- function(arms, n, a.time = NULL, a.rate = NULL, d.hazard = 0,
                           endpoint = c("os", "pfs"), alpha = 0.025, power = 0.8,
                           tmax = NULL) {
  endpoint <- match.arg(endpoint)
  if (!is.list(arms) || inherits(arms, "OncoArm") || length(arms) != 2L) {
    stop("'arms' must be a list of two OncoArm objects (control first).", call. = FALSE)
  }
  if (!(alpha > 0 && alpha < 0.5) || !(power > 0 && power < 1)) {
    stop("'alpha' must be in (0, 0.5) and 'power' in (0, 1).", call. = FALSE)
  }
  if (is.null(tmax)) {
    a_end <- if (is.null(a.time)) 0 else max(a.time)
    meds <- vapply(.check_arms(arms), QuantileEndpoint, numeric(1), prob = 0.5,
                   endpoint = endpoint)
    tmax <- a_end + 10 * max(meds)
  }
  setup <- .design_setup(arms, n, a.time, a.rate, d.hazard, endpoint, tmax = tmax)
  ev <- function(tcal) .design_integrals(setup, tcal)$events
  e_max <- ev(tmax)
  r1 <- n[1] / sum(n)
  z2 <- (stats::qnorm(1 - alpha) + stats::qnorm(power)) ^ 2
  tcal <- stats::uniroot(function(x) ev(x) - 0.5 * e_max, c(1e-8, tmax), tol = 1e-10)$root
  a <- .design_integrals(setup, tcal)$ahr
  d_prev <- NA_integer_
  it <- 0L
  repeat {
    it <- it + 1L
    if (abs(log(a)) < 1e-12) stop("The average hazard ratio is 1.", call. = FALSE)
    d_exact <- z2 / (r1 * (1 - r1) * log(a) ^ 2)
    d <- as.integer(ceiling(d_exact))
    if (d > e_max) {
      stop("The required number of events (", d, ") exceeds the expected ",
           "number of events by 'tmax' (", round(e_max, 1), ").", call. = FALSE)
    }
    tcal <- stats::uniroot(function(x) ev(x) - d, c(1e-8, tmax), tol = 1e-10)$root
    a_new <- .design_integrals(setup, tcal)$ahr
    converged <- identical(d, d_prev) && abs(a_new - a) < 1e-10
    a <- a_new
    if (converged || it >= 100L) break
    d_prev <- d
  }
  if (!converged) warning("The iteration did not converge in 100 steps.", call. = FALSE)
  list(events = d, events.exact = d_exact, time = tcal, ahr = a,
       max.events = e_max, iterations = it)
}
