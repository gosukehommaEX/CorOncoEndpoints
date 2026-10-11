#' Check and tabulate a piecewise-uniform accrual distribution (internal)
#'
#' Accrual is piecewise uniform over the intervals defined by \code{a.time}.
#' The expected numbers of patients in the intervals are given by absolute
#' rates (\code{a.rate}, patients per unit time over all groups, as in the
#' FastSurvival package), by proportions (\code{a.prop}) or, when neither is
#' given, by the lengths of the intervals. With absolute rates, either the
#' rates accrue exactly \code{n_total} patients over the intervals, or one
#' more rate than intervals is given and the last rate applies to an open
#' final interval whose end is computed from \code{n_total}.
#'
#' @param a.time Accrual breakpoints (non-negative and increasing), or
#'   \code{NULL} when all patients enter at time 0.
#' @param a.rate Absolute accrual rates, \code{length(a.time) - 1} or
#'   \code{length(a.time)} values.
#' @param a.prop Proportions of patients in the intervals,
#'   \code{length(a.time) - 1} values.
#' @param n_total Total number of patients over all groups.
#' @return A list with elements \code{time} (breakpoints, including a computed
#'   end), \code{prob} (probabilities of the intervals) and \code{cum}
#'   (cumulative probabilities at the breakpoints).
#' @keywords internal
#' @noRd
.accrual_setup <- function(a.time, a.rate = NULL, a.prop = NULL, n_total) {
  if (is.null(a.time)) {
    if (!is.null(a.rate) || !is.null(a.prop)) {
      stop("'a.rate' and 'a.prop' require 'a.time'.", call. = FALSE)
    }
    return(list(time = 0, prob = numeric(0), cum = 0))
  }
  if (!is.numeric(a.time) || length(a.time) < 1L || any(!is.finite(a.time)) ||
      a.time[1] < 0 || any(diff(a.time) <= 0)) {
    stop("'a.time' must be a strictly increasing finite numeric vector of ",
         "non-negative times.", call. = FALSE)
  }
  if (!is.null(a.rate) && !is.null(a.prop)) {
    stop("Give at most one of 'a.rate' and 'a.prop'.", call. = FALSE)
  }
  k <- length(a.time) - 1L
  if (!is.null(a.rate)) {
    if (!is.numeric(a.rate) || any(!is.finite(a.rate)) || any(a.rate < 0)) {
      stop("'a.rate' must contain non-negative accrual rates.", call. = FALSE)
    }
    if (length(a.rate) == k + 1L) {
      if (!(a.rate[k + 1L] > 0)) {
        stop("The rate of the open final accrual interval must be positive.",
             call. = FALSE)
      }
      bounded <- if (k > 0L) sum(a.rate[seq_len(k)] * diff(a.time)) else 0
      remaining <- n_total - bounded
      if (!(remaining > 0)) {
        stop("The accrual intervals of 'a.time' already accrue ", signif(bounded, 6),
             " patients, at least the total of 'n' (", n_total, "); the open ",
             "final interval cannot be added.", call. = FALSE)
      }
      a.time <- c(a.time, a.time[k + 1L] + remaining / a.rate[k + 1L])
      w <- a.rate * diff(a.time)
    } else if (length(a.rate) == k && k > 0L) {
      w <- a.rate * diff(a.time)
      if (abs(sum(w) - n_total) > 1e-8 * max(1, n_total)) {
        stop("'a.rate' accrues ", signif(sum(w), 6), " patients over 'a.time' ",
             "but the total of 'n' is ", n_total, ". Make them consistent, give ",
             "one more rate for an open final interval, or use 'a.prop' for ",
             "relative accrual.", call. = FALSE)
      }
    } else {
      stop("'a.rate' must have length(a.time) - 1 elements (accrual over the ",
           "intervals of 'a.time') or length(a.time) elements (the last rate ",
           "applies to an open final interval).", call. = FALSE)
    }
  } else {
    if (k < 1L) {
      stop("'a.time' must have at least two elements unless 'a.rate' gives the ",
           "rate of an open final interval.", call. = FALSE)
    }
    if (is.null(a.prop)) {
      w <- diff(a.time)
    } else {
      if (!is.numeric(a.prop) || length(a.prop) != k || any(!is.finite(a.prop)) ||
          any(a.prop < 0) || sum(a.prop) <= 0) {
        stop("'a.prop' must have length(a.time) - 1 non-negative elements with a ",
             "positive sum.", call. = FALSE)
      }
      w <- as.numeric(a.prop)
    }
  }
  if (!(sum(w) > 0)) stop("The accrual rates must have a positive sum.", call. = FALSE)
  list(time = a.time, prob = w / sum(w), cum = c(0, cumsum(w) / sum(w)))
}
