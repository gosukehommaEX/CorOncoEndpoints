#' Generate correlated PFS, OS and response
#'
#' @description
#' Generates patient-level progression-free survival (PFS), overall survival
#' (OS), objective response and, optionally, time to response for one or more
#' treatment groups and many simulated trials, with accrual and dropout.
#'
#' @details
#' Each group is described by an \code{\link{OncoArm}} object; see its Details
#' for the model. Random numbers come from the \code{dqrng} generator and the
#' transformation runs in C++. For reproducible results set \code{seed}; the
#' \code{dqrng} generator is not affected by \code{set.seed()}.
#'
#' Accrual is piecewise uniform: patients enter independently, with density
#' proportional to \code{a.rate} on the intervals defined by \code{a.time}.
#' Without \code{a.time} all patients enter at time 0. Dropout is exponential
#' with hazard \code{d.hazard} and censors PFS, OS and the observation of
#' response.
#'
#' The latent times \code{pfs_time} and \code{os_time} always satisfy
#' \code{pfs_time <= os_time}, with equality when the PFS event is a death
#' (\code{progression = 0}). The observed columns \code{pfs_tte},
#' \code{pfs_event}, \code{os_tte} and \code{os_event} account for dropout but
#' not for an analysis cutoff; use \code{\link{CutoffData}} for that.
#'
#' @param nsim Number of simulated trials.
#' @param n Integer vector of sample sizes, one per group.
#' @param arms An \code{OncoArm} object or a list of them, one per group. Group
#'   labels are taken from the names of the list, otherwise from the
#'   \code{label} of each object, otherwise \code{"Group1"}, \code{"Group2"},
#'   and so on.
#' @param a.time Accrual breakpoints starting at 0, or \code{NULL}.
#' @param a.rate Relative accrual intensities on the intervals of
#'   \code{a.time} (default: equal).
#' @param d.hazard Dropout hazard, a scalar or one value per group.
#' @param seed Optional integer seed for \code{dqrng::dqset.seed()}.
#' @return A data frame with one row per patient, ordered by \code{sim} and
#'   group, with columns \code{sim}, \code{group}, \code{accrual_time},
#'   \code{pfs_time}, \code{os_time}, \code{progression} (1 if the PFS event is
#'   a progression, 0 if a death), \code{response}, \code{ttr} (time to
#'   response, \code{NA} unless \code{resp.timing = "ttr"} and the patient
#'   responded), \code{dropout_time}, \code{pfs_tte}, \code{pfs_event},
#'   \code{os_tte}, \code{os_event}, \code{pfs_calendar_time} and
#'   \code{os_calendar_time}.
#' @examples
#' ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
#'                death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6)
#' trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
#'                death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
#' dat <- rOncoEndpoints(nsim = 2, n = c(100, 100),
#'                       arms = list(Control = ctl, Treatment = trt),
#'                       a.time = c(0, 24), seed = 1)
#' head(dat)
#' @seealso \code{\link{OncoArm}}, \code{\link{CutoffData}}
#' @export
rOncoEndpoints <- function(nsim = 1, n, arms, a.time = NULL, a.rate = NULL,
                           d.hazard = 0, seed = NULL) {
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim < 1 || nsim != round(nsim)) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }
  arms <- .check_arms(arms)
  k <- length(arms)
  if (!is.numeric(n) || length(n) != k || any(n < 1) || any(n != round(n))) {
    stop("'n' must be a vector of positive integers with one element per group.",
         call. = FALSE)
  }
  if (!is.numeric(d.hazard) || !(length(d.hazard) %in% c(1L, k)) ||
      any(!is.finite(d.hazard)) || any(d.hazard < 0)) {
    stop("'d.hazard' must be a non-negative number or one per group.", call. = FALSE)
  }
  if (nsim * max(n) > .Machine$integer.max) {
    stop("'nsim' times the largest group size must not exceed ",
         .Machine$integer.max, "; generate the trials in batches.", call. = FALSE)
  }
  d.hazard <- rep_len(d.hazard, k)
  acc <- .accrual_setup(a.time, a.rate)
  if (!is.null(seed)) dqrng::dqset.seed(seed)

  parts <- vector("list", k)
  for (j in seq_len(k)) {
    a <- arms[[j]]
    expexp <- a$os.model == "expexp"
    out <- generate_arm_cpp(
      nsim = as.integer(nsim), n = as.integer(n[j]),
      lam_p = a$lam_p, theta = a$theta, c_resp = a$c_resp, z_tau = a$z_tau,
      tau = a$tau, os_model = if (expexp) 1L else 0L,
      pi_death = a$death.prop,
      gam0 = if (expexp) 1 else a$gam0, gam1 = if (expexp) 1 else a$gam1,
      lam_o = if (expexp) a$lam_o else 1, c_dec = if (expexp) a$c_dec else 0,
      grid_g = if (expexp) a$grid$G else c(0, 0),
      grid_h = if (expexp) a$grid$h else 1,
      ttr_mode = if (a$resp.timing == "ttr") 1L else 0L,
      ttr_rate = if (a$resp.timing == "ttr") a$ttr_rate else 1,
      a_time = acc$time, a_cum = acc$cum, d_hazard = d.hazard[j])
    out$group <- rep(names(arms)[j], length(out$sim))
    out$group_index <- rep(j, length(out$sim))
    parts[[j]] <- out
  }
  col <- function(name) unlist(lapply(parts, `[[`, name), use.names = FALSE)
  res <- data.frame(sim = col("sim"), group = col("group"),
                    accrual_time = col("accrual_time"),
                    pfs_time = col("pfs_time"), os_time = col("os_time"),
                    progression = col("progression"), response = col("response"),
                    ttr = col("ttr"), dropout_time = col("dropout_time"),
                    stringsAsFactors = FALSE)
  ord <- order(res$sim, col("group_index"))
  res <- res[ord, , drop = FALSE]
  res$pfs_tte <- pmin(res$pfs_time, res$dropout_time)
  res$pfs_event <- as.integer(res$pfs_time <= res$dropout_time)
  res$os_tte <- pmin(res$os_time, res$dropout_time)
  res$os_event <- as.integer(res$os_time <= res$dropout_time)
  res$pfs_calendar_time <- res$accrual_time + res$pfs_tte
  res$os_calendar_time <- res$accrual_time + res$os_tte
  rownames(res) <- NULL
  res
}
