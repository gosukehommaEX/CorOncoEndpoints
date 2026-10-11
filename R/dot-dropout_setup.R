#' Check and tabulate the dropout distribution (internal)
#'
#' Dropout is exponential, or piecewise exponential when breakpoints are given
#' in \code{d.time}. Without \code{d.time}, \code{d.hazard} (or
#' \code{d.median}) is a scalar or has one value per group. With
#' \code{d.time}, a numeric \code{d.hazard} gives the hazards of the pieces for
#' all groups and a list gives them per group; \code{d.time} is a numeric
#' vector for all groups or a list with one element per group, starting at 0
#' and ending with \code{Inf}.
#'
#' @param d.hazard,d.median,d.time See \code{rOncoEndpoints()}.
#' @param k Number of groups.
#' @return A list with one element per group, each a list with elements
#'   \code{time} (left ends of the pieces, starting at 0), \code{hazard} and
#'   \code{cum} (cumulative hazard at \code{time}).
#' @keywords internal
#' @noRd
.dropout_setup <- function(d.hazard = 0, d.median = NULL, d.time = NULL, k) {
  if (!is.null(d.median)) {
    no_hazard <- is.numeric(d.hazard) && length(d.hazard) == 1L &&
      isTRUE(d.hazard == 0)
    if (!no_hazard) {
      stop("Give at most one of 'd.hazard' and 'd.median'.", call. = FALSE)
    }
    to_hazard <- function(m) {
      if (!is.numeric(m) || length(m) < 1L || anyNA(m) || any(m <= 0)) {
        stop("'d.median' must contain positive medians (Inf for no dropout).",
             call. = FALSE)
      }
      log(2) / m
    }
    d.hazard <- if (is.list(d.median)) lapply(d.median, to_hazard) else to_hazard(d.median)
  }
  if (is.list(d.hazard)) {
    if (length(d.hazard) != k) {
      stop("A list 'd.hazard' or 'd.median' must have one element per group.",
           call. = FALSE)
    }
    hz <- d.hazard
  } else if (is.null(d.time)) {
    if (!is.numeric(d.hazard) || !(length(d.hazard) %in% c(1L, k))) {
      stop("'d.hazard' must be a non-negative number or one per group.", call. = FALSE)
    }
    hz <- as.list(rep_len(d.hazard, k))
  } else {
    hz <- rep(list(d.hazard), k)
  }
  if (is.list(d.time)) {
    if (length(d.time) != k) {
      stop("A list 'd.time' must have one element per group.", call. = FALSE)
    }
    tm <- d.time
  } else {
    tm <- rep(list(d.time), k)
  }
  lapply(seq_len(k), function(j) {
    h <- hz[[j]]
    if (!is.numeric(h) || length(h) < 1L || any(!is.finite(h)) || any(h < 0)) {
      stop("'d.hazard' must be a non-negative number or one per group.", call. = FALSE)
    }
    b <- tm[[j]]
    if (is.null(b)) {
      if (length(h) != 1L) {
        stop("Piecewise dropout hazards require 'd.time'.", call. = FALSE)
      }
      return(list(time = 0, hazard = as.numeric(h), cum = 0))
    }
    if (!is.numeric(b) || length(b) < 2L || anyNA(b) || b[1] != 0 ||
        b[length(b)] != Inf || any(diff(b) <= 0)) {
      stop("'d.time' must be increasing, start at 0 and end with Inf.", call. = FALSE)
    }
    if (length(h) != length(b) - 1L) {
      stop("'d.hazard' must have length(d.time) - 1 elements.", call. = FALSE)
    }
    left <- b[-length(b)]
    list(time = left, hazard = as.numeric(h),
         cum = c(0, cumsum(h[-length(h)] * diff(left))))
  })
}
