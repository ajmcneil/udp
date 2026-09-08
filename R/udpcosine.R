# udp transformations based on cosine functions.

#' Class of cosine udp transformations
#'
#' A cosine udp transformation is the uniform-distribution-preserving map of
#' `[0, 1]` given by a folded cosine of the stated degree: for `U` uniform on
#' `[0, 1]`, `udptrans(x, U)` is again uniform on `[0, 1]`.
#'
#' @slot degree integer; the degree of the transformation (at least 1).
#'
#' @seealso [udpcosine()] to construct one, [udptrans()] to evaluate it.
#' @include udp-package.R
#' @export
setClass("udpcosine", contains = "udp", slots = list(degree = "integer"))

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

#' @describeIn udptrans Evaluate a cosine udp transformation.
#' @export
setMethod("udptrans", "udpcosine", function(x, u) {
  d <- x@degree
  arg <- pmin(pmax(cos(d * pi * u) * (-1)^d, -1), 1)
  out <- 1 - acos(arg) / pi
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' Roots of a cosine udp transformation
#'
#' `udptrans()` is a degree-fold triangle wave and so is not injective: a
#' value `v` in `(0, 1)` has `degree` pre-images, one in each linear piece,
#' while `v = 0` and `v = 1` have fewer (the shared troughs and peaks
#' respectively).
#'
#' @param x an object of class \linkS4class{udpcosine}.
#' @param v a vector with values in `[0, 1]`.
#'
#' @return A list the same length as `v`; element `j` is the sorted vector of
#'   `u` in `[0, 1]` with `udptrans(x, u)` equal to `v[j]`.
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

#' @describeIn udpsi Stochastic inverse of a cosine udp transformation. For each
#'   `v` in `(0, 1)`, [udpcosinverse()] returns `degree` pre-images and the
#'   selected one is root number `ceiling(Z * degree)`; `v` equal to `0` or `1`
#'   returns `0`.
#' @export
setMethod("udpsi", "udpcosine", function(x, v, Z = runif(length(v)), ...) {
  if (anyNA(v) || any(v < 0 | v > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  if (length(Z) != length(v)) {
    stop("'Z' must have the same length as 'v'.", call. = FALSE)
  }
  d <- x@degree
  edge <- v <= 0 | v >= 1
  vw <- v
  vw[edge] <- 0.5 # working value with the full complement of degree roots
  j <- pmin(pmax(ceiling(as.numeric(Z) * d), 1L), d)

  if (length(v) == 0L) {
    out <- numeric(0)
  } else {
    roots <- do.call(rbind, udpcosinverse(x, vw)) # length(v) x degree
    out <- roots[cbind(seq_along(v), j)]
  }
  out[edge] <- 0

  if (!is.null(attributes(v))) {
    attributes(out) <- attributes(v)
  }
  out
})

#' Plot method for the udpcosine class
#'
#' Draws the graph of the cosine udp transformation over thin red gridlines at
#' the boundaries of its `degree` linear pieces.
#'
#' @param x an object of class \linkS4class{udpcosine}.
#' @param xlab,ylab axis labels.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(udpcosine(3))
#' plot(udpcosine(4))
setMethod("plot", c(x = "udpcosine", y = "missing"),
  function(x, xlab = "u", ylab = "T(u)", ...) {
    u <- (0:x@degree) / x@degree
    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
      xlab = xlab, ylab = ylab, ...
    )
    abline(v = u, col = "red", lwd = 0.5)
    lines(u, udptrans(x, u), lwd = 2)
  }
)
