#' Gauss-Legendre nodes and weights (internal)
#'
#' Computes the nodes and weights of the \code{n}-point Gauss-Legendre rule on
#' \code{[-1, 1]} from the eigen decomposition of the Jacobi matrix (Golub and
#' Welsch, 1969).
#'
#' @param n Number of nodes.
#' @return A list with elements \code{x} (increasing nodes) and \code{w}
#'   (weights).
#' @keywords internal
#' @noRd
.gauss_legendre <- function(n) {
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k ^ 2 - 1)
  jac <- matrix(0, n, n)
  jac[cbind(k, k + 1L)] <- b
  jac[cbind(k + 1L, k)] <- b
  e <- eigen(jac, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = 2 * e$vectors[1, ord] ^ 2)
}
