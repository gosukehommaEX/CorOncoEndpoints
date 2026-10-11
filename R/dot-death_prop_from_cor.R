#' Proportion of deaths that gives a target correlation of PFS and OS (internal)
#'
#' Finds \code{death.prop} such that the Pearson correlation between PFS and OS
#' of the group equals \code{target}, by root finding over \code{range}. The
#' correlation decreases as \code{death.prop} increases.
#'
#' @param target Target \code{Corr(PFS, OS)}.
#' @param build Function of \code{death.prop} that returns the \code{OncoArm}
#'   object with all other inputs fixed.
#' @param range Search interval of \code{death.prop}.
#' @return Scalar \code{death.prop}.
#' @keywords internal
#' @noRd
.death_prop_from_cor <- function(target, build, range) {
  rho <- function(dp) CorEndpoints(build(dp))[["cor.pfs.os"]]
  r_lo <- rho(range[1])
  r_hi <- rho(range[2])
  if (target > max(r_lo, r_hi) || target < min(r_lo, r_hi)) {
    stop("'pfs.os.cor' = ", target, " is outside the attainable range (",
         round(min(r_lo, r_hi), 4), ", ", round(max(r_lo, r_hi), 4),
         ") for these inputs, obtained with 'death.prop' from ",
         signif(range[1], 3), " to ", signif(range[2], 3), ".", call. = FALSE)
  }
  stats::uniroot(function(dp) rho(dp) - target, interval = range,
                 f.lower = r_lo - target, f.upper = r_hi - target,
                 tol = 1e-10)$root
}
