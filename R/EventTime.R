#' Calendar time of a given number of events
#'
#' Returns, for each simulated trial, the calendar time at which the given
#' number of observed PFS or OS events is reached. Events censored by dropout
#' are not counted. The result can be passed to \code{\link{CutoffData}} for
#' event-driven analyses.
#'
#' @param data A data frame returned by \code{\link{rOncoEndpoints}}.
#' @param events Number of events.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A numeric vector with one value per simulated trial (\code{NA} when
#'   fewer events occur), named by \code{sim}.
#' @examples
#' arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15)
#' dat <- rOncoEndpoints(nsim = 3, n = 200, arms = arm, a.time = c(0, 12),
#'                       seed = 1)
#' EventTime(dat, events = 100, endpoint = "os")
#' @export
EventTime <- function(data, events, endpoint = c("os", "pfs")) {
  endpoint <- match.arg(endpoint)
  if (!is.numeric(events) || length(events) != 1L || events < 1 ||
      events != round(events)) {
    stop("'events' must be a positive integer.", call. = FALSE)
  }
  ev_col <- paste0(endpoint, "_event")
  cal_col <- paste0(endpoint, "_calendar_time")
  if (!is.data.frame(data) || !all(c("sim", ev_col, cal_col) %in% names(data))) {
    stop("'data' must be a data frame returned by rOncoEndpoints().", call. = FALSE)
  }
  ev <- data[[ev_col]] == 1
  sims <- sort(unique(data$sim))
  sp <- split(data[[cal_col]][ev], factor(data$sim[ev], levels = sims))
  vapply(sp, function(x) {
    if (length(x) >= events) sort(x, partial = events)[events] else NA_real_
  }, numeric(1))
}
