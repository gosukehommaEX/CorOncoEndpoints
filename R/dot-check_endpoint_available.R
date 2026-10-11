#' Check that a group has the requested endpoint (internal)
#'
#' @param arm An \code{OncoArm} object.
#' @param endpoint \code{"pfs"} or \code{"os"}.
#' @param response \code{"all"}, \code{"responders"} or \code{"nonresponders"};
#'   the last two require response among the endpoints.
#' @return \code{TRUE}, invisibly, or an error.
#' @keywords internal
#' @noRd
.check_endpoint_available <- function(arm, endpoint, response = "all") {
  ep <- .arm_endpoints(arm)
  if (!endpoint %in% ep) {
    stop("The group has no ", toupper(endpoint), " (endpoints: ",
         paste(ep, collapse = ", "), ").", call. = FALSE)
  }
  if (response != "all" && !"response" %in% ep) {
    stop("The group has no response, so 'response' must be \"all\".", call. = FALSE)
  }
  invisible(TRUE)
}
