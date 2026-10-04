#' Expected events and average hazard ratio at a calendar time (internal)
#'
#' Uses the midpoint rule on the intervals of the tabulation grid up to the
#' calendar time \code{tcal}. The weight of follow-up time \code{x} is the
#' expected number of events at \code{x} among patients who entered at least
#' \code{x} before \code{tcal}:
#' \code{n_j f_j(x) exp(-d_j x) A(tcal - x)}, where \code{A} is the accrual
#' distribution function.
#'
#' @param setup Output of \code{.design_setup()}.
#' @param tcal Calendar time.
#' @return A list with the expected events per group (\code{by_group}), their
#'   sum (\code{events}) and, for two groups, the average hazard ratio of the
#'   second group versus the first (\code{ahr}).
#' @keywords internal
#' @noRd
.design_integrals <- function(setup, tcal) {
  g1 <- setup$grids[[1]]
  tt <- g1$t
  if (tcal > tt[length(tt)]) stop("Calendar time beyond the tabulation grid.", call. = FALSE)
  edges <- c(tt[tt < tcal], tcal)
  if (length(edges) < 2L) edges <- c(0, tcal)
  mid <- 0.5 * (edges[-1] + edges[-length(edges)])
  width <- diff(edges)
  a_w <- .accrual_cdf(tcal - mid, setup$acc)
  k <- length(setup$grids)
  w_list <- vector("list", k)
  h_list <- vector("list", k)
  by_group <- numeric(k)
  for (j in seq_len(k)) {
    g <- setup$grids[[j]]
    fm <- stats::approx(g$t, g$f, xout = mid)$y
    hm <- stats::approx(g$t, g$h, xout = mid)$y
    w_list[[j]] <- setup$n[j] * fm * exp(-setup$d[j] * mid) * a_w
    h_list[[j]] <- hm
    by_group[j] <- sum(w_list[[j]] * width)
  }
  names(by_group) <- names(setup$arms)
  ahr <- NA_real_
  if (k == 2L) {
    w <- (w_list[[1]] + w_list[[2]]) * width
    lhr <- log(h_list[[2]] / h_list[[1]])
    ok <- is.finite(lhr) & w > 0
    ahr <- exp(sum(lhr[ok] * w[ok]) / sum(w[ok]))
  }
  list(by_group = by_group, events = sum(by_group), ahr = ahr)
}
