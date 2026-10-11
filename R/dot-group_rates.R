#' Check a non-negative rate given for all groups or per group (internal)
#'
#' @param x Numeric scalar or vector with one value per group.
#' @param k Number of groups.
#' @param name Name of the argument, for error messages.
#' @return Numeric vector of length \code{k}.
#' @keywords internal
#' @noRd
.group_rates <- function(x, k, name) {
  if (!is.numeric(x) || !(length(x) %in% c(1L, k)) || any(!is.finite(x)) ||
      any(x < 0)) {
    stop("'", name, "' must be a non-negative number or one per group.", call. = FALSE)
  }
  rep_len(as.numeric(x), k)
}
