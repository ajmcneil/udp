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

#' Roots of a cosine udp transformation
#'
#' `udpcostrans()` is a degree-fold triangle wave and so is not injective: a
#' value `v` in `(0, 1)` has `degree` pre-images, one in each linear piece,
#' while `v = 0` and `v = 1` have fewer (the shared troughs and peaks
#' respectively).
#'
#' @param x an object of class \linkS4class{udpcosine}.
#' @param v a vector with values in `[0, 1]`.
#'
#' @return A list the same length as `v`; element `j` is the sorted vector of
#'   `u` in `[0, 1]` with `udpcostrans(x, u)` equal to `v[j]`.
#' @export
#'
#' @examples
#' udpcosinverse(udpcosine(3), c(0, 0.4, 1))
udpcosinverse <- function(x, v) {
  if (anyNA(v) || any(v < 0 | v > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  d <- x@degree
  k <- seq_len(d)
  increasing <- (k %% 2L) == (d %% 2L)
  lapply(v, function(vi) {
    sort(unique(ifelse(increasing, (vi + k - 1) / d, (k - vi) / d)))
  })
}
