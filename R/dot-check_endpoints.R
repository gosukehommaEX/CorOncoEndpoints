#' Check the endpoints of a group (internal)
#'
#' @param endpoints Character vector with two or three of \code{"pfs"},
#'   \code{"os"} and \code{"response"}.
#' @return The endpoints in the order \code{"pfs"}, \code{"os"},
#'   \code{"response"}, without duplicates.
#' @keywords internal
#' @noRd
.check_endpoints <- function(endpoints) {
  all3 <- c("pfs", "os", "response")
  if (!is.character(endpoints) || anyNA(endpoints) || !all(endpoints %in% all3)) {
    stop("'endpoints' must contain \"pfs\", \"os\" or \"response\".", call. = FALSE)
  }
  ep <- all3[all3 %in% endpoints]
  if (length(ep) < 2L) {
    stop("'endpoints' must contain at least two of \"pfs\", \"os\" and \"response\".",
         call. = FALSE)
  }
  ep
}
