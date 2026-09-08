# udp transformations based on cosine functions.

#' Class of cosine udp transformations
#'
#' A cosine udp transformation is the uniform-distribution-preserving map of
#' `[0, 1]` given by a folded cosine of the stated degree: for `U` uniform on
#' `[0, 1]`, `udpcostrans(x, U)` is again uniform on `[0, 1]`.
#'
#' @slot degree integer; the degree of the transformation (at least 1).
#'
#' @seealso [udpcosine()] to construct one, [udpcostrans()] to evaluate it.
#' @export
setClass("udpcosine", slots = list(degree = "integer"))

#' Construct a cosine udp transformation
#'
#' @param degree a single positive integer, the degree of the transformation.
#'
#' @return An object of class \linkS4class{udpcosine}.
#' @export
#'
#' @examples
#' udpcosine(2)
udpcosine <- function(degree) {
  if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) ||
    degree < 1 || degree != round(degree)) {
    stop("'degree' must be a single positive integer.", call. = FALSE)
  }
  new("udpcosine", degree = as.integer(degree))
}

#' Evaluate a cosine udp transformation
#'
#' @param x an object of class \linkS4class{udpcosine}.
#' @param u a vector with values in `[0, 1]`.
#'
#' @return An object shaped like `u` with values in `[0, 1]`.
#' @export
#'
#' @examples
#' udpcostrans(udpcosine(2), c(0, 0.25, 0.5, 0.75, 1))
udpcostrans <- function(x, u) {
  d <- x@degree
  arg <- pmin(pmax(cos(d * pi * u) * (-1)^d, -1), 1)
  out <- 1 - acos(arg) / pi
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
}
