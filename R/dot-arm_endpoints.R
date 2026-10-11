#' Endpoints of a group (internal)
#'
#' @param arm An \code{OncoArm} object; objects created by version 0.2.0 have
#'   no \code{endpoints} element and have all three endpoints.
#' @return Character vector of endpoints.
#' @keywords internal
#' @noRd
.arm_endpoints <- function(arm) {
  if (is.null(arm$endpoints)) c("pfs", "os", "response") else arm$endpoints
}
