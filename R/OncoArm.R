#' Define a treatment group for correlated PFS, OS and response
#'
#' @description
#' Creates an \code{OncoArm} object that describes one treatment group and
#' calibrates the parameters of the generator to design inputs: the PFS median
#' (or hazard), the objective response rate, the association between PFS and
#' response, the proportion of PFS events that are deaths, and the OS median or
#' the post-progression survival.
#'
#' @details
#' \strong{PFS and response.} Let \code{(Z1, Z2)} be standard bivariate normal
#' with correlation \code{theta}. PFS is \code{-log(Phi(-Z1)) / lam_p}, so it
#' is exactly exponential with hazard \code{lam_p}. Response is
#' \code{R = 1} when \code{PFS > resp.tau} and \code{Z2 > c_resp}, where
#' \code{c_resp} is chosen so that \code{P(R = 1) = orr} exactly. The copula
#' correlation \code{theta} is calibrated either to the Pearson correlation
#' \code{resp.cor} between PFS and response or to the median PFS among
#' responders \code{resp.pfs.median}. Exactly one of these two must be given.
#'
#' \strong{Timing of response.} With \code{resp.timing = "none"} response is a
#' patient-level binary variable without a time. With
#' \code{resp.timing = "landmark"} responders must have PFS longer than
#' \code{resp.tau} (for example, the time of the first tumor assessment). With
#' \code{resp.timing = "ttr"} the same restriction applies and, in addition,
#' each responder gets a time to response that lies between \code{resp.tau}
#' and the responder's PFS: \code{resp.tau} plus an exponential time with
#' median \code{ttr.median - resp.tau}, truncated at \code{PFS - resp.tau}.
#' The time to response allows the response observed by an interim analysis
#' to be derived (see \code{\link{CutoffData}}). In all three modes, PFS
#' remains exactly exponential and the response probability remains exactly
#' \code{orr}.
#'
#' \strong{OS, illness-death model (\code{os.model = "idm"}).} A PFS event is a
#' death with probability \code{death.prop}, independently of PFS and response;
#' the pre-progression hazards of progression and death are therefore
#' constant, \code{(1 - death.prop) * lam_p} and \code{death.prop * lam_p}.
#' After progression, the patient survives for an exponential time with hazard
#' \code{gam0} (non-responders) or \code{gam1 = pps.hr.resp * gam0}
#' (responders), measured from progression. \code{gam0} is calibrated to
#' \code{os.median} or given through \code{pps.median} or \code{pps.hazard}.
#' OS is not exponential unless \code{pps.hr.resp = 1} and
#' \code{gam0 = death.prop * lam_p}.
#'
#' \strong{OS, exponential-exponential model (\code{os.model = "expexp"}).} OS
#' is exactly exponential with median \code{os.median}. This requires the
#' pre-progression death hazard to start at the OS hazard \code{lam_o}; it is
#' set to \code{lam_o * exp(-c_dec * t)} with \code{c_dec} chosen so that the
#' proportion of PFS events that are deaths equals \code{death.prop}, which
#' must not exceed \code{pfs.median / os.median}. The post-progression hazard
#' is then determined by the exponential OS distribution (it depends on time
#' since randomization) and cannot depend on response, so
#' \code{pps.hr.resp} must be 1. With
#' \code{death.prop = pfs.median / os.median} the model reduces to the maximal
#' independence model of Fleischer et al. (2009).
#'
#' @param pfs.median,pfs.hazard PFS median or hazard; give exactly one.
#' @param orr Objective response rate, in (0, 1).
#' @param resp.cor Target Pearson correlation between PFS and response.
#' @param resp.pfs.median Target median PFS among responders. Give exactly one
#'   of \code{resp.cor} and \code{resp.pfs.median}.
#' @param death.prop Proportion of PFS events that are deaths without prior
#'   progression, in [0, 1) for \code{"idm"} and in (0, pfs.median /
#'   os.median] for \code{"expexp"}.
#' @param os.model \code{"idm"} (default) or \code{"expexp"}; see Details.
#' @param os.median Target OS median. Required for \code{"expexp"}; for
#'   \code{"idm"} give exactly one of \code{os.median}, \code{pps.median} and
#'   \code{pps.hazard}.
#' @param pps.median,pps.hazard Median or hazard of post-progression survival
#'   of non-responders (\code{"idm"} only).
#' @param pps.hr.resp Post-progression hazard ratio of responders versus
#'   non-responders (\code{"idm"} only; must be 1 for \code{"expexp"}).
#' @param resp.timing \code{"none"} (default), \code{"landmark"} or
#'   \code{"ttr"}; see Details.
#' @param resp.tau Landmark time: responders have PFS longer than
#'   \code{resp.tau}. Must be 0 for \code{"none"} and positive for
#'   \code{"landmark"}; for \code{"ttr"} it may be 0 (no landmark).
#' @param ttr.median Median of the untruncated time to response, larger than
#'   \code{resp.tau} (\code{"ttr"} only).
#' @param label Optional group label used by \code{\link{rOncoEndpoints}}.
#'
#' @return An object of class \code{OncoArm}: a list with the inputs and the
#'   calibrated parameters \code{lam_p} (PFS hazard), \code{theta} (copula
#'   correlation), \code{c_resp} (response threshold), \code{z_tau} (latent
#'   score of the landmark), and for \code{"idm"} \code{gam0} and \code{gam1},
#'   or for \code{"expexp"} \code{lam_o}, \code{c_dec} and a tabulated
#'   cumulative post-progression hazard.
#'
#' @references
#' Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A statistical
#' model for the dependence between progression-free survival and overall
#' survival. \emph{Statistics in Medicine}, 28, 2669--2686.
#' \doi{10.1002/sim.3637}
#'
#' @examples
#' ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
#'                death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6,
#'                label = "Control")
#' ctl
#'
#' # Both PFS and OS exactly exponential
#' ee <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
#'               death.prop = 0.15, os.model = "expexp", os.median = 15)
#' QuantileEndpoint(ee, endpoint = "os")
#'
#' @seealso \code{\link{rOncoEndpoints}}, \code{\link{CorEndpoints}}
#' @export
OncoArm <- function(pfs.median = NULL, pfs.hazard = NULL, orr,
                    resp.cor = NULL, resp.pfs.median = NULL,
                    death.prop,
                    os.model = c("idm", "expexp"),
                    os.median = NULL, pps.median = NULL, pps.hazard = NULL,
                    pps.hr.resp = 1,
                    resp.timing = c("none", "landmark", "ttr"),
                    resp.tau = 0, ttr.median = NULL,
                    label = NULL) {
  os.model <- match.arg(os.model)
  resp.timing <- match.arg(resp.timing)
  is_pos <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0

  # ---- PFS ----
  if (is.null(pfs.median) == is.null(pfs.hazard)) {
    stop("Give exactly one of 'pfs.median' and 'pfs.hazard'.", call. = FALSE)
  }
  if (!is.null(pfs.median)) {
    if (!is_pos(pfs.median)) stop("'pfs.median' must be a positive number.", call. = FALSE)
    lam_p <- log(2) / pfs.median
  } else {
    if (!is_pos(pfs.hazard)) stop("'pfs.hazard' must be a positive number.", call. = FALSE)
    lam_p <- pfs.hazard
  }
  pfs_med <- log(2) / lam_p

  # ---- response ----
  if (!is.numeric(orr) || length(orr) != 1L || is.na(orr) || !(orr > 0 && orr < 1)) {
    stop("'orr' must be a number in (0, 1).", call. = FALSE)
  }
  if (is.null(resp.cor) == is.null(resp.pfs.median)) {
    stop("Give exactly one of 'resp.cor' and 'resp.pfs.median'.", call. = FALSE)
  }
  if (!is.numeric(resp.tau) || length(resp.tau) != 1L || !is.finite(resp.tau) ||
      resp.tau < 0) {
    stop("'resp.tau' must be a non-negative number.", call. = FALSE)
  }
  if (resp.timing == "none" && resp.tau != 0) {
    stop("'resp.tau' must be 0 when resp.timing = \"none\".", call. = FALSE)
  }
  if (resp.timing == "landmark" && resp.tau <= 0) {
    stop("'resp.tau' must be positive when resp.timing = \"landmark\".", call. = FALSE)
  }
  if (resp.timing == "ttr") {
    if (!is_pos(ttr.median) || ttr.median <= resp.tau) {
      stop("'ttr.median' must be a number larger than 'resp.tau' when ",
           "resp.timing = \"ttr\".", call. = FALSE)
    }
  } else if (!is.null(ttr.median)) {
    stop("'ttr.median' is used only when resp.timing = \"ttr\".", call. = FALSE)
  }
  if (resp.tau > 0 && orr >= exp(-lam_p * resp.tau)) {
    stop("'orr' must be smaller than the probability that PFS exceeds ",
         "'resp.tau' (", round(exp(-lam_p * resp.tau), 4), ").", call. = FALSE)
  }
  z_tau <- if (resp.tau > 0) .z_of_time(resp.tau, lam_p) else -Inf
  if (!is.null(resp.cor)) {
    if (!is.numeric(resp.cor) || length(resp.cor) != 1L || is.na(resp.cor) ||
        abs(resp.cor) >= 1) {
      stop("'resp.cor' must be a number in (-1, 1).", call. = FALSE)
    }
    theta <- .theta_from_cor(resp.cor, orr, z_tau)
  } else {
    if (!is_pos(resp.pfs.median) || resp.pfs.median <= resp.tau) {
      stop("'resp.pfs.median' must be a number larger than 'resp.tau'.", call. = FALSE)
    }
    theta <- .theta_from_resp_median(resp.pfs.median, orr, lam_p, z_tau)
  }
  c_resp <- .resp_threshold(theta, orr, z_tau)

  # ---- OS ----
  if (!is.numeric(death.prop) || length(death.prop) != 1L || is.na(death.prop) ||
      death.prop < 0 || death.prop >= 1) {
    stop("'death.prop' must be a number in [0, 1).", call. = FALSE)
  }
  if (!is_pos(pps.hr.resp)) stop("'pps.hr.resp' must be a positive number.", call. = FALSE)
  if (!is.null(os.median)) {
    if (!is_pos(os.median)) stop("'os.median' must be a positive number.", call. = FALSE)
    if (os.median <= pfs_med) {
      stop("'os.median' must be larger than the PFS median.", call. = FALSE)
    }
  }

  arm <- list(os.model = os.model, resp.timing = resp.timing, label = label,
              lam_p = lam_p, orr = orr, theta = theta, c_resp = c_resp,
              z_tau = z_tau, tau = resp.tau, death.prop = death.prop,
              pps.hr.resp = pps.hr.resp,
              ttr_rate = if (resp.timing == "ttr") log(2) / (ttr.median - resp.tau) else NA_real_,
              inputs = list(pfs.median = pfs.median, pfs.hazard = pfs.hazard,
                            resp.cor = resp.cor, resp.pfs.median = resp.pfs.median,
                            os.median = os.median, pps.median = pps.median,
                            pps.hazard = pps.hazard, ttr.median = ttr.median))
  class(arm) <- "OncoArm"

  if (os.model == "idm") {
    n_given <- sum(!is.null(os.median), !is.null(pps.median), !is.null(pps.hazard))
    if (n_given != 1L) {
      stop("For os.model = \"idm\" give exactly one of 'os.median', ",
           "'pps.median' and 'pps.hazard'.", call. = FALSE)
    }
    if (!is.null(os.median)) {
      gam0 <- .calib_gam0(arm, os.median)
    } else if (!is.null(pps.median)) {
      if (!is_pos(pps.median)) stop("'pps.median' must be a positive number.", call. = FALSE)
      gam0 <- log(2) / pps.median
    } else {
      if (!is_pos(pps.hazard)) stop("'pps.hazard' must be a positive number.", call. = FALSE)
      gam0 <- pps.hazard
    }
    arm$gam0 <- gam0
    arm$gam1 <- pps.hr.resp * gam0
  } else {
    if (is.null(os.median)) {
      stop("'os.median' is required for os.model = \"expexp\".", call. = FALSE)
    }
    if (!is.null(pps.median) || !is.null(pps.hazard)) {
      stop("'pps.median' and 'pps.hazard' cannot be used with ",
           "os.model = \"expexp\".", call. = FALSE)
    }
    if (pps.hr.resp != 1) {
      stop("'pps.hr.resp' must be 1 for os.model = \"expexp\".", call. = FALSE)
    }
    lam_o <- log(2) / os.median
    if (!(death.prop > 0) || death.prop > lam_o / lam_p * (1 + 1e-12)) {
      stop("For os.model = \"expexp\", 'death.prop' must be in (0, ",
           round(lam_o / lam_p, 6), "], the ratio of the PFS median to the ",
           "OS median.", call. = FALSE)
    }
    arm$lam_o <- lam_o
    arm$c_dec <- max(lam_o / death.prop - lam_p, 0)
    arm$grid <- .expexp_grid(lam_p, lam_o, arm$c_dec)
  }
  arm
}
