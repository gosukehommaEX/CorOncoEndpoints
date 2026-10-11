#' Generate correlated PFS, OS and response
#'
#' @description
#' Generates patient-level progression-free survival (PFS), overall survival
#' (OS), objective response and, optionally, time to response for one or more
#' treatment groups and many simulated trials, with accrual, dropout and
#' censoring of PFS.
#'
#' @details
#' Each group is described by an \code{\link{OncoArm}} object; see its Details
#' for the model. Random numbers come from the \code{dqrng} generator and the
#' transformation runs in C++. For reproducible results set \code{seed}; the
#' \code{dqrng} generator is not affected by \code{set.seed()}.
#'
#' \strong{Accrual.} Accrual is piecewise uniform over the intervals defined
#' by \code{a.time}, as in the FastSurvival package. Each group has a fixed
#' number of patients in each interval: the group size is divided among the
#' intervals in proportion to their expected numbers of patients, with the
#' remainder going to the intervals with the largest fractional parts, and the
#' patients are placed uniformly within their interval. The expected numbers
#' are given in one of three ways.
#' \itemize{
#'   \item \code{a.rate}: absolute rates, in patients per unit time over all
#'     groups. With \code{length(a.time) - 1} rates, the rates must accrue the
#'     total of \code{n}. With \code{length(a.time)} rates, the last rate applies
#'     to an open final interval whose end is computed so that the total of
#'     \code{n} is reached; for example, \code{a.time = 0} with
#'     \code{a.rate = 20} enrolls 20 patients per unit time until the total is
#'     reached.
#'   \item \code{a.prop}: the proportions of patients in the intervals.
#'   \item Neither: proportional to the lengths of the intervals, that is,
#'     uniform accrual from the first to the last element of \code{a.time}.
#' }
#' \code{a.time} may start after 0, for example for patients enrolled in the
#' expansion of a trial. Without \code{a.time} all patients enter at time 0.
#'
#' \strong{Dropout and censoring of PFS.} Dropout is exponential with hazard
#' \code{d.hazard} (or median \code{d.median}), or piecewise exponential with
#' breakpoints \code{d.time}, and censors PFS and OS. \code{pfs.cens.hazard}
#' adds an exponential censoring time that censors PFS only, for example the
#' start of another anticancer therapy before progression, while OS remains
#' under follow-up. Both are independent of the endpoints. Whether a response
#' is observed before dropout, before the censoring of PFS and before an
#' analysis cutoff is determined by \code{\link{CutoffData}}.
#'
#' \strong{Random numbers.} The groups draw their random numbers in turn. The
#' censoring time of PFS is drawn only for groups with a positive
#' \code{pfs.cens.hazard}, so that data generated without it do not change.
#' With \code{stream}, the generator is seeded with \code{seed} and moved to an
#' independent stream, so that batches of simulated trials can be generated
#' with the same \code{seed} and different streams, in sequence or in
#' parallel.
#'
#' The latent times \code{pfs_time} and \code{os_time} always satisfy
#' \code{pfs_time <= os_time}, with equality when the PFS event is a death
#' (\code{progression = 0}). The observed columns \code{pfs_tte},
#' \code{pfs_event}, \code{os_tte} and \code{os_event} account for dropout and
#' the censoring of PFS but not for an analysis cutoff, and \code{response} is
#' the response of the patient whether or not it would be observed; use
#' \code{\link{CutoffData}} for the data observed at an analysis. For groups
#' with two endpoints (\code{endpoints} in \code{\link{OncoArm}}), the columns
#' of the missing endpoint are \code{NA}; all groups must have the same
#' endpoints.
#'
#' @param nsim Number of simulated trials.
#' @param n Integer vector of sample sizes, one per group.
#' @param arms An \code{OncoArm} object or a list of them, one per group. Group
#'   labels are taken from the names of the list, otherwise from the
#'   \code{label} of each object, otherwise \code{"Group1"}, \code{"Group2"},
#'   and so on.
#' @param a.time Accrual breakpoints (non-negative and increasing), or
#'   \code{NULL} for entry of all patients at time 0.
#' @param a.rate Absolute accrual rates (patients per unit time over all
#'   groups), with \code{length(a.time) - 1} or \code{length(a.time)}
#'   elements; see Details.
#' @param a.prop Proportions of patients enrolled in the intervals of
#'   \code{a.time}, with \code{length(a.time) - 1} elements. Give at most one
#'   of \code{a.rate} and \code{a.prop}.
#' @param d.hazard Dropout hazard: a scalar or one value per group, or, with
#'   \code{d.time}, the hazards of the pieces (a numeric vector for all groups
#'   or a list with one vector per group).
#' @param d.median Dropout median, an alternative to \code{d.hazard} with the
#'   same structure (\code{Inf} for no dropout).
#' @param d.time Breakpoints of a piecewise exponential dropout distribution,
#'   starting at 0 and ending with \code{Inf}: a numeric vector for all groups
#'   or a list with one vector per group.
#' @param pfs.cens.hazard Hazard of an exponential censoring time that censors
#'   PFS only, a scalar or one value per group (0 for none).
#' @param seed Optional integer seed for \code{dqrng::dqset.seed()}.
#' @param stream Optional non-negative whole number that selects an independent
#'   stream of the generator for \code{seed}; requires \code{seed}.
#' @return A data frame with one row per patient, ordered by \code{sim} and
#'   group, with columns \code{sim}, \code{group}, \code{accrual_time},
#'   \code{pfs_time}, \code{os_time}, \code{progression} (1 if the PFS event is
#'   a progression, 0 if a death), \code{response}, \code{ttr} (time to
#'   response, \code{NA} unless \code{resp.timing = "ttr"} and the patient
#'   responded), \code{dropout_time}, \code{pfs_tte}, \code{pfs_event},
#'   \code{os_tte}, \code{os_event}, \code{pfs_calendar_time},
#'   \code{os_calendar_time}, \code{pfs_cens_time} (censoring time of PFS,
#'   \code{Inf} without \code{pfs.cens.hazard}) and \code{resp_tau} (landmark
#'   of the group, used by \code{\link{CutoffData}}).
#' @examples
#' ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
#'                death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6)
#' trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
#'                death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
#' dat <- rOncoEndpoints(nsim = 2, n = c(100, 100),
#'                       arms = list(Control = ctl, Treatment = trt),
#'                       a.time = c(0, 24), seed = 1)
#' head(dat)
#'
#' # 10 patients per month for 6 months, then 30 per month until 200 patients
#' # are enrolled (month 6 + 140 / 30); dropout of 5% per year (months)
#' dat2 <- rOncoEndpoints(nsim = 2, n = c(100, 100),
#'                        arms = list(Control = ctl, Treatment = trt),
#'                        a.time = c(0, 6), a.rate = c(10, 30),
#'                        d.hazard = -log(0.95) / 12, seed = 1)
#' range(dat2$accrual_time)
#' @seealso \code{\link{OncoArm}}, \code{\link{CutoffData}}
#' @export
rOncoEndpoints <- function(nsim = 1, n, arms, a.time = NULL, a.rate = NULL,
                           a.prop = NULL, d.hazard = 0, d.median = NULL,
                           d.time = NULL, pfs.cens.hazard = 0, seed = NULL,
                           stream = NULL) {
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim < 1 || nsim != round(nsim)) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }
  arms <- .check_arms(arms)
  k <- length(arms)
  if (!is.numeric(n) || length(n) != k || any(n < 1) || any(n != round(n))) {
    stop("'n' must be a vector of positive integers with one element per group.",
         call. = FALSE)
  }
  eps <- lapply(arms, .arm_endpoints)
  if (!all(vapply(eps, identical, logical(1), eps[[1]]))) {
    stop("All groups must have the same 'endpoints'.", call. = FALSE)
  }
  ep <- eps[[1]]
  # 0: all three endpoints; 1: no response; 2: no OS; 3: no PFS
  ep_mode <- if (length(ep) == 3L) 0L else if (!"response" %in% ep) 1L else
    if (!"os" %in% ep) 2L else 3L
  if (nsim * max(n) > .Machine$integer.max) {
    stop("'nsim' times the largest group size must not exceed ",
         .Machine$integer.max, "; generate the trials in batches.", call. = FALSE)
  }
  acc <- .accrual_setup(a.time, a.rate, a.prop, sum(n))
  dsp <- .dropout_setup(d.hazard, d.median, d.time, k)
  pc <- .group_rates(pfs.cens.hazard, k, "pfs.cens.hazard")
  if (!"pfs" %in% ep && any(pc > 0)) {
    stop("'pfs.cens.hazard' requires PFS in 'endpoints'.", call. = FALSE)
  }
  if (!is.null(stream)) {
    if (is.null(seed)) stop("'stream' requires 'seed'.", call. = FALSE)
    if (!is.numeric(stream) || length(stream) != 1L || !is.finite(stream) ||
        stream < 0 || stream != round(stream) || stream > .Machine$integer.max) {
      stop("'stream' must be a non-negative whole number.", call. = FALSE)
    }
    stream <- as.integer(stream)
  }
  if (!is.null(seed)) dqrng::dqset.seed(seed, stream = stream)

  parts <- vector("list", k)
  for (j in seq_len(k)) {
    a <- arms[[j]]
    idm <- a$os.model == "idm"
    expexp <- a$os.model == "expexp"
    out <- generate_arm_cpp(
      nsim = as.integer(nsim), n = as.integer(n[j]),
      lam_p = a$lam_p, theta = a$theta, c_resp = a$c_resp, z_tau = a$z_tau,
      tau = a$tau, os_model = if (expexp) 1L else 0L,
      pi_death = if (idm || expexp) a$death.prop else 0,
      gam0 = if (idm) a$gam0 else 1, gam1 = if (idm) a$gam1 else 1,
      lam_o = if (expexp) a$lam_o else 1, c_dec = if (expexp) a$c_dec else 0,
      grid_g = if (expexp) a$grid$G else c(0, 0),
      grid_h = if (expexp) a$grid$h else 1,
      ttr_mode = if (a$resp.timing == "ttr") 1L else 0L,
      ttr_rate = if (a$resp.timing == "ttr") a$ttr_rate else 1,
      a_time = acc$time, a_counts = .accrual_counts(n[j], acc$prob),
      d_time = dsp[[j]]$time, d_hazard = dsp[[j]]$hazard, d_cum = dsp[[j]]$cum,
      pc_hazard = pc[j], ep_mode = ep_mode)
    parts[[j]] <- out
  }
  # Rows are ordered by trial and, within a trial, by group. Each row of group j
  # is written directly to its position instead of sorting all rows.
  n_all <- sum(n)
  first <- c(0, cumsum(n))[seq_len(k)]
  pos <- lapply(seq_len(k), function(j) {
    rep((seq_len(nsim) - 1) * n_all, each = n[j]) + rep.int(first[j] + seq_len(n[j]), nsim)
  })
  col <- function(name) {
    x <- parts[[1L]][[name]]
    if (k == 1L) return(x)
    out <- vector(typeof(x), nsim * n_all)
    for (j in seq_len(k)) out[pos[[j]]] <- parts[[j]][[name]]
    out
  }
  res <- data.frame(sim = rep(seq_len(nsim), each = n_all),
                    group = rep(rep(names(arms), n), nsim),
                    accrual_time = col("accrual_time"),
                    pfs_time = col("pfs_time"), os_time = col("os_time"),
                    progression = col("progression"), response = col("response"),
                    ttr = col("ttr"), dropout_time = col("dropout_time"),
                    stringsAsFactors = FALSE)
  pfs_cens <- col("pfs_cens_time")
  pfs_fu <- if (any(pc > 0)) pmin(res$dropout_time, pfs_cens) else res$dropout_time
  res$pfs_tte <- pmin(res$pfs_time, pfs_fu)
  res$pfs_event <- as.integer(res$pfs_time <= pfs_fu)
  res$os_tte <- pmin(res$os_time, res$dropout_time)
  res$os_event <- as.integer(res$os_time <= res$dropout_time)
  res$pfs_calendar_time <- res$accrual_time + res$pfs_tte
  res$os_calendar_time <- res$accrual_time + res$os_tte
  res$pfs_cens_time <- pfs_cens
  taus <- vapply(arms, function(a) as.numeric(a$tau), numeric(1), USE.NAMES = FALSE)
  res$resp_tau <- rep(rep(unname(taus), n), nsim)
  rownames(res) <- NULL
  res
}
