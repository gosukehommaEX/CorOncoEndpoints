#' Survival, density and hazard of PFS or OS
#'
#' @description
#' Evaluates the survival function, density or hazard of PFS or OS implied by an
#' \code{OncoArm} object, for all patients or separately for responders and
#' non-responders.
#'
#' @details
#' PFS is exponential for all patients; among responders its survival function
#' is a bivariate normal orthant probability divided by the response rate. For
#' OS with \code{os.model = "idm"} and \code{pps.hr.resp = 1}, the survival
#' function of all patients has the closed form of Fleischer et al. (2009,
#' Theorem 5); otherwise the contribution of patients who progressed is a
#' one-dimensional integral over the progression time. For
#' \code{os.model = "expexp"}, OS of all patients is exponential.
#'
#' @param arm An object of class \code{OncoArm}.
#' @param t Numeric vector of non-negative times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param response \code{"all"}, \code{"responders"} or \code{"nonresponders"}.
#' @param type \code{"survival"}, \code{"density"} or \code{"hazard"}.
#' @return Numeric vector of the same length as \code{t}.
#' @references
#' Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A statistical
#' model for the dependence between progression-free survival and overall
#' survival. \emph{Statistics in Medicine}, 28, 2669--2686.
#' \doi{10.1002/sim.3637}
#' @examples
#' arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'                os.median = 15, pps.hr.resp = 0.6)
#' SurvEndpoint(arm, c(6, 12, 24), endpoint = "os")
#' SurvEndpoint(arm, c(6, 12, 24), endpoint = "os", response = "responders")
#' SurvEndpoint(arm, c(6, 12, 24), endpoint = "os", type = "hazard")
#' @export
SurvEndpoint <- function(arm, t, endpoint = c("os", "pfs"),
                         response = c("all", "responders", "nonresponders"),
                         type = c("survival", "density", "hazard")) {
  if (!inherits(arm, "OncoArm")) stop("'arm' must be an OncoArm object.", call. = FALSE)
  endpoint <- match.arg(endpoint)
  response <- match.arg(response)
  type <- match.arg(type)
  if (!is.numeric(t) || any(!is.finite(t)) || any(t < 0)) {
    stop("'t' must be a numeric vector of non-negative finite times.", call. = FALSE)
  }
  if (endpoint == "pfs") {
    return(.surv_pfs(t, arm, response = response, type = type))
  }
  if (type == "hazard") {
    return(.surv_os(t, arm, response, "density") / .surv_os(t, arm, response, "survival"))
  }
  .surv_os(t, arm, response, type)
}
