#' Compute Correlation Bounds Between Response and Time-to-Event Endpoints
#'
#' @description
#' Computes the lower and upper bounds of the correlation coefficient between
#' a binary response endpoint and a time-to-event (TTE) endpoint. These bounds
#' are derived using Fréchet-Hoeffding copula bounds and depend only on the
#' response probability, not on the hazard rate of the TTE endpoint.
#'
#' @param p Numeric. The true probability of the binary response endpoint.
#'   Must be between 0 and 1 (exclusive).
#'
#' @return A numeric vector of length 2 containing:
#'   \itemize{
#'     \item First element: Lower bound (maximum negative dependence)
#'     \item Second element: Upper bound (maximum positive dependence)
#'   }
#'
#' @details
#' The Fréchet-Hoeffding bounds provide theoretical limits for the correlation
#' between a binary response endpoint (e.g., objective response) and an
#' exponentially distributed TTE endpoint (e.g., overall survival or
#' progression-free survival).
#'
#' The correlation bounds are given by:
#' \deqn{\log(1-p)\sqrt{\frac{1-p}{p}} \leq \text{Corr}(TTE, Response) \leq \log\left(\frac{1}{p}\right)\sqrt{\frac{p}{1-p}}}
#'
#' Key properties:
#' \itemize{
#'   \item The bounds depend only on \code{p}, not on the TTE hazard rate
#'   \item Lower bound: achieved by countermonotonic copula \eqn{C(u,v) = \max\{u+v-1, 0\}}
#'   \item Upper bound: achieved by comonotonic copula \eqn{C(u,v) = \min\{u, v\}}
#'   \item The bounds are symmetric: if \code{p} is replaced by \code{1-p},
#'         the lower and upper bounds swap with negated signs
#'   \item Maximum absolute bound value of approximately 0.805 occurs at
#'         \code{p ≈ 0.203} (for upper bound) and \code{p ≈ 0.797} (for lower bound).
#'         This is a fundamental property of the Fréchet-Hoeffding bounds for
#'         correlations between binary and continuous exponential variables.
#' }
#'
#' @note
#' The independence of correlation bounds from the hazard rate is an
#' interesting theoretical result. While the hazard rate determines the scale
#' of the TTE endpoint, it does not affect the range of achievable correlations
#' with the binary response endpoint. This occurs because the hazard rate only
#' affects the scaling of the TTE variable, while the binary response endpoint
#' is determined solely by the probability \code{p}.
#'
#' @references
#' Trivedi, P. K., & Zimmer, D. M. (2005). Copula modeling: an introduction
#' for practitioners. Foundations and Trends in Econometrics, 1(1), 1-111.
#'
#' @examples
#' # Calculate bounds with response probability = 0.4
#' CorBoundResponseTTE(p = 0.4)
#'
#' # Example with higher response rate
#' CorBoundResponseTTE(p = 0.6)
#'
#' # Example with lower response rate
#' CorBoundResponseTTE(p = 0.2)
#'
#' # Example with intermediate values
#' bounds <- CorBoundResponseTTE(p = 0.35)
#' cat("Lower bound:", bounds[1], "\nUpper bound:", bounds[2], "\n")
#'
#' # Demonstrating the effect of response probability on bounds
#' # The bounds are widest around p ≈ 0.203 and p ≈ 0.797
#' CorBoundResponseTTE(p = 0.203)
#' CorBoundResponseTTE(p = 0.5)
#' CorBoundResponseTTE(p = 0.797)
#'
#' # Demonstrating symmetry property
#' CorBoundResponseTTE(p = 0.3)
#' -rev(CorBoundResponseTTE(p = 0.7))  # Should be nearly identical
#'
#' # Find maximum absolute bound
#' p_seq <- seq(0.01, 0.99, by = 0.001)
#' bounds_upper <- sapply(p_seq, function(p) CorBoundResponseTTE(p)[2])
#' max_bound_idx <- which.max(bounds_upper)
#' cat("Maximum upper bound:", round(bounds_upper[max_bound_idx], 4),
#'     "at p =", round(p_seq[max_bound_idx], 3), "\n")
#'
#' @export
CorBoundResponseTTE <- function(p) {
  # Input validation
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("'p' must be a single numeric value between 0 and 1 (exclusive)")
  }
  
  # Calculate correlation bounds using Fréchet-Hoeffding copula bounds
  # Lower bound: log(1-p) * sqrt((1-p)/p) using countermonotonic copula
  # Upper bound: log(1/p) * sqrt(p/(1-p)) using comonotonic copula
  # Note: The bounds are independent of the TTE hazard rate
  rho_lower <- log(1 - p) * sqrt((1 - p) / p)
  rho_upper <- log(1 / p) * sqrt(p / (1 - p))
  
  rho_bound <- c(rho_lower, rho_upper)
  
  return(rho_bound)
}