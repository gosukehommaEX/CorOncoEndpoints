#' Quantiles of PFS or OS
#'
#' Computes the time at which the survival function of PFS or OS, for all
#' patients or by response, equals \code{1 - prob}. The median is
#' \code{prob = 0.5}.
#'
#' @param arm An object of class \code{OncoArm}.
#' @param prob Numeric vector of probabilities in (0, 1).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param response \code{"all"}, \code{"responders"} or \code{"nonresponders"}.
#' @return Numeric vector of times.
#' @examples
#' arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15)
#' QuantileEndpoint(arm, endpoint = "os")
#' QuantileEndpoint(arm, endpoint = "pfs", response = "responders")
#' @export
QuantileEndpoint <- function(arm, prob = 0.5, endpoint = c("os", "pfs"),
                             response = c("all", "responders", "nonresponders")) {
  if (!inherits(arm, "OncoArm")) stop("'arm' must be an OncoArm object.", call. = FALSE)
  endpoint <- match.arg(endpoint)
  response <- match.arg(response)
  if (!is.numeric(prob) || anyNA(prob) || any(!(prob > 0 & prob < 1))) {
    stop("'prob' must be in (0, 1).", call. = FALSE)
  }
  vapply(prob, function(pr) {
    if (endpoint == "pfs" && response == "all") return(-log(1 - pr) / arm$lam_p)
    if (endpoint == "os" && response == "all" && arm$os.model == "expexp") {
      return(-log(1 - pr) / arm$lam_o)
    }
    S <- function(tt) SurvEndpoint(arm, tt, endpoint = endpoint, response = response)
    target <- 1 - pr
    hi <- max(-log(1 - pr) / arm$lam_p, arm$tau) * 2 + 1
    while (S(hi) > target) hi <- hi * 2
    stats::uniroot(function(tt) S(tt) - target, interval = c(0, hi), tol = 1e-10)$root
  }, numeric(1))
}
