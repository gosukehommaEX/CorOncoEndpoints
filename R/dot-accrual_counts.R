#' Numbers of patients per accrual interval (internal)
#'
#' Distributes \code{n} patients over the accrual intervals in proportion to
#' \code{prob}: each interval receives the integer part of its share, and the
#' remaining patients go to the intervals with the largest fractional parts
#' (largest-remainder method, as in the FastSurvival package), so that the
#' counts add up to \code{n}. Ties go to the earlier interval.
#'
#' @param n Number of patients (a whole number).
#' @param prob Probabilities of the intervals, summing to one; empty when all
#'   patients enter at the same time.
#' @return Integer vector of the same length as \code{prob}.
#' @keywords internal
#' @noRd
.accrual_counts <- function(n, prob) {
  if (length(prob) == 0L) return(integer(0))
  target <- n * prob
  base <- floor(target)
  rem <- as.integer(round(n - sum(base)))
  if (rem > 0L) {
    take <- order(target - base, decreasing = TRUE)[seq_len(rem)]
    base[take] <- base[take] + 1
  }
  as.integer(base)
}
