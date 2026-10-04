#' Check a list of OncoArm objects and their sample sizes (internal)
#'
#' @param arms A single \code{OncoArm} object or a list of them.
#' @return A named list of \code{OncoArm} objects.
#' @keywords internal
#' @noRd
.check_arms <- function(arms) {
  if (inherits(arms, "OncoArm")) arms <- list(arms)
  if (!is.list(arms) || length(arms) < 1L ||
      !all(vapply(arms, inherits, logical(1), what = "OncoArm"))) {
    stop("'arms' must be an OncoArm object or a list of OncoArm objects.",
         call. = FALSE)
  }
  nm <- names(arms)
  if (is.null(nm)) nm <- rep("", length(arms))
  for (j in seq_along(arms)) {
    if (nm[j] == "") {
      nm[j] <- if (!is.null(arms[[j]]$label)) arms[[j]]$label else paste0("Group", j)
    }
  }
  if (anyDuplicated(nm)) stop("Group labels must be unique.", call. = FALSE)
  names(arms) <- nm
  arms
}
