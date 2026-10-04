#' Attainable range of the correlation between PFS and response
#'
#' @description
#' Returns the smallest and largest Pearson correlation between an exponential
#' PFS and a binary response with probability \code{orr}, optionally under the
#' restriction that responders have PFS longer than a landmark time
#' \code{resp.tau}.
#'
#' @details
#' Without a landmark the bounds are the Frechet bounds
#' \code{sqrt((1 - p) / p) log(1 - p)} and \code{-sqrt(p / (1 - p)) log(p)},
#' which do not depend on the PFS hazard. With a landmark the upper bound is
#' unchanged, and the lower bound is attained when the responders are the
#' patients with the shortest PFS above \code{resp.tau}. The Gaussian copula
#' used by \code{\link{OncoArm}} approaches both bounds as the copula
#' correlation tends to -1 or 1; \code{OncoArm()} restricts the copula
#' correlation to (-0.9999, 0.9999), so correlations very close to a bound
#' may not be attainable.
#'
#' @param orr Objective response rate, in (0, 1).
#' @param pfs.median,pfs.hazard PFS median or hazard; needed only when
#'   \code{resp.tau > 0}.
#' @param resp.tau Landmark time (0 for none).
#' @return A numeric vector \code{c(lower = , upper = )}.
#' @examples
#' CorBoundPFSResponse(orr = 0.3)
#' CorBoundPFSResponse(orr = 0.3, pfs.median = 6, resp.tau = 1.5)
#' @export
CorBoundPFSResponse <- function(orr, pfs.median = NULL, pfs.hazard = NULL,
                                resp.tau = 0) {
  if (!is.numeric(orr) || length(orr) != 1L || !(orr > 0 && orr < 1)) {
    stop("'orr' must be a number in (0, 1).", call. = FALSE)
  }
  p <- orr
  upper <- -sqrt(p / (1 - p)) * log(p)
  if (resp.tau == 0) {
    lower <- sqrt((1 - p) / p) * log(1 - p)
  } else {
    if (is.null(pfs.median) == is.null(pfs.hazard)) {
      stop("Give exactly one of 'pfs.median' and 'pfs.hazard' when ",
           "'resp.tau' > 0.", call. = FALSE)
    }
    lam <- if (!is.null(pfs.hazard)) pfs.hazard else log(2) / pfs.median
    a <- lam * resp.tau
    eb <- exp(-a) - p
    if (eb <= 0) {
      stop("'orr' must be smaller than the probability that PFS exceeds ",
           "'resp.tau'.", call. = FALSE)
    }
    lb <- -log(eb)
    lower <- ((a + 1) * exp(-a) - (lb + 1) * eb - p) / sqrt(p * (1 - p))
  }
  c(lower = lower, upper = upper)
}
