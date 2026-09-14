# GENERAL-PURPOSE UTILITIES NOT TIED TO ANY PARTICULAR udp CLASS

#' Force Values Strictly Inside (0, 1)
#'
#' Nudges any values exactly at the boundary of `[0, 1]` strictly inside,
#' replacing `0` by `btol` and `1` by `1 - btol`. Useful whenever a fitted
#' model or transformation produces values that should lie in the open
#' interval `(0, 1)` but can attain the closed boundary exactly -- for
#' example the output of [udptrans()], or of a marginal distribution fitted
#' to data and then evaluated at the data themselves.
#'
#' @param u a numeric vector, matrix or array with values in `[0, 1]`.
#' @param btol boundary tolerance. If `NA` (the default), `1 / (2 * n)` is
#'   used, where `n` is the number of rows of `u` if `u` is a matrix (one
#'   tolerance per observation, however many columns it has) and
#'   `length(u)` otherwise.
#'
#' @return `u` with any values exactly at `0` or `1` replaced so that all
#'   values lie strictly inside `(0, 1)`.
#' @export
#'
#' @examples
#' boundaryadjust(c(0, 0.3, 1))
#' boundaryadjust(matrix(c(0, 0.4, 1, 0.7), 2, 2))
boundaryadjust <- function(u, btol = NA) {
  if (is.na(btol)) {
    n <- if (is.matrix(u)) nrow(u) else length(u)
    btol <- 1 / (2 * n)
  }
  u[u == 0] <- btol
  u[u == 1] <- 1 - btol
  u
}
