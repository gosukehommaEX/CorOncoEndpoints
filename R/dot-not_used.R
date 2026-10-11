#' Error for arguments that the endpoints of a group do not use (internal)
#'
#' The arguments \code{cond} and \code{why} follow \code{...}, so that they
#' are matched only by their full names and an argument to check cannot be
#' taken for them by partial matching.
#'
#' @param ... Named arguments to check.
#' @param cond Logical; when \code{TRUE}, the arguments in \code{...} must be
#'   \code{NULL}.
#' @param why Text that completes the error message.
#' @return \code{NULL}, invisibly, or an error naming the first argument given.
#' @keywords internal
#' @noRd
.not_used <- function(..., cond, why) {
  if (!cond) return(invisible(NULL))
  vals <- list(...)
  given <- names(vals)[!vapply(vals, is.null, logical(1))]
  if (length(given) > 0L) {
    stop("'", given[1], "' is not used ", why, ".", call. = FALSE)
  }
  invisible(NULL)
}
