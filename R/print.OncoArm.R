#' Print an OncoArm object
#'
#' Prints the calibrated parameters of a treatment group together with the
#' implied medians of PFS and OS (overall and by response) and the implied
#' Pearson correlations among PFS, OS and response.
#'
#' @param x An object of class \code{OncoArm}.
#' @param digits Number of significant digits.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @examples
#' OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
#'         os.median = 15)
#' @export
print.OncoArm <- function(x, digits = 4, ...) {
  f <- function(v) format(signif(v, digits))
  lab <- if (is.null(x$label)) "" else paste0(" (", x$label, ")")
  cat("OncoArm", lab, "\n", sep = "")
  cat("  OS model           : ", x$os.model, "\n", sep = "")
  cat("  Response timing    : ", x$resp.timing,
      if (x$tau > 0) paste0(" (tau = ", f(x$tau), ")") else "", "\n", sep = "")
  cat("  PFS hazard         : ", f(x$lam_p), "\n", sep = "")
  cat("  Response rate      : ", f(x$orr), "\n", sep = "")
  cat("  Copula correlation : ", f(x$theta), "\n", sep = "")
  cat("  Death proportion   : ", f(x$death.prop), "\n", sep = "")
  if (x$os.model == "idm") {
    cat("  Post-progression hazard (non-responders, responders): ",
        f(x$gam0), ", ", f(x$gam1), "\n", sep = "")
  } else {
    cat("  OS hazard          : ", f(x$lam_o), "\n", sep = "")
    cat("  Decay of pre-progression death hazard: ", f(x$c_dec), "\n", sep = "")
  }
  med <- function(endpoint, response) {
    QuantileEndpoint(x, 0.5, endpoint = endpoint, response = response)
  }
  cat("Implied medians (all, responders, non-responders)\n")
  cat("  PFS: ", f(med("pfs", "all")), ", ", f(med("pfs", "responders")), ", ",
      f(med("pfs", "nonresponders")), "\n", sep = "")
  cat("  OS : ", f(med("os", "all")), ", ", f(med("os", "responders")), ", ",
      f(med("os", "nonresponders")), "\n", sep = "")
  cr <- CorEndpoints(x)
  cat("Implied correlations\n")
  cat("  Corr(PFS, R) = ", f(cr[["cor.pfs.resp"]]),
      ", Corr(OS, R) = ", f(cr[["cor.os.resp"]]),
      ", Corr(PFS, OS) = ", f(cr[["cor.pfs.os"]]), "\n", sep = "")
  cat("  Partial correlation of OS and R given PFS = ",
      f(cr[["pcor.os.resp"]]), "\n", sep = "")
  invisible(x)
}
