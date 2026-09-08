# Shuffles: uniform-distribution-preserving, piecewise-linear bijections of
# [0, 1] with slopes of 1 and -1.

#' Class of shuffles
#'
#' A shuffle splits `[0, 1]` into `m` equal strips and maps domain strip `i`
#' onto range strip `perm[i]`: order-preserving (slope `1`) where `signs[i]` is
#' `1`, order-reversing (slope `-1`) where it is `-1`. The result is a
#' uniform-distribution-preserving bijection of `[0, 1]`.
#'
#' @slot perm integer vector; a permutation of `seq_len(m)`.
#' @slot signs numeric vector of `1` and `-1`, the same length as `perm`.
#'
#' @seealso [shuffle()] to construct one; [shtrans()] and [shinverse()] to
#'   evaluate it and its inverse.
#' @export
setClass("shuffle", slots = list(perm = "integer", signs = "numeric"))

#' Construct a shuffle
#'
#' @param perm a permutation of `1:m`, given as a numeric or integer vector.
#' @param signs a vector of `1` and `-1` of length `m` giving the slope of the
#'   shuffle on each domain strip. Defaults to all `1` (no flips).
#'
#' @return An object of class \linkS4class{shuffle}.
#' @export
#'
#' @examples
#' shuffle(c(2, 1, 3))
#' shuffle(c(3, 1, 2), signs = c(1, -1, 1))
shuffle <- function(perm, signs = rep(1, length(perm))) {
  perm <- as.integer(perm)
  m <- length(perm)
  if (m < 1L || anyNA(perm) || !setequal(perm, seq_len(m))) {
    stop("'perm' must be a permutation of 1:m.", call. = FALSE)
  }
  if (length(signs) != m || anyNA(signs) || !all(signs %in% c(-1, 1))) {
    stop(
      "'signs' must be a vector of 1 and -1 the same length as 'perm'.",
      call. = FALSE
    )
  }
  new("shuffle", perm = perm, signs = as.numeric(signs))
}

#' Evaluate a shuffle
#'
#' @param x an object of class \linkS4class{shuffle}.
#' @param u a vector with values in `[0, 1]`.
#'
#' @return An object shaped like `u` with values in `[0, 1]`.
#' @export
#'
#' @examples
#' s <- shuffle(c(2, 1, 3))
#' shtrans(s, c(0, 0.2, 0.5, 0.9, 1))
shtrans <- function(x, u) {
  m <- length(x@perm)
  i <- as.integer(pmax(pmin(floor(u * m) + 1, m), 1))
  t <- u * m - (i - 1)
  s <- x@signs[i]
  out <- (x@perm[i] - (s > 0) + s * t) / m
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
}

#' Evaluate the inverse of a shuffle
#'
#' The inverse of a shuffle is itself a shuffle, with permutation
#' `order(perm)` and the signs reordered to match.
#'
#' @param x an object of class \linkS4class{shuffle}.
#' @param v a vector with values in `[0, 1]`.
#'
#' @return An object shaped like `v` with values in `[0, 1]`.
#' @export
#'
#' @examples
#' s <- shuffle(c(3, 1, 2), signs = c(1, -1, 1))
#' shinverse(s, shtrans(s, c(0.1, 0.5, 0.9)))
shinverse <- function(x, v) {
  o <- order(x@perm)
  shtrans(new("shuffle", perm = o, signs = x@signs[o]), v)
}

#' Plot method for the shuffle class
#'
#' Draws the graph of the shuffle as one thick black line segment per strip,
#' over thin red gridlines at the boundaries of the vertical and horizontal
#' strips.
#'
#' @param x an object of class \linkS4class{shuffle}.
#' @param xlab,ylab axis labels.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(shuffle(c(3, 1, 2), signs = c(1, -1, 1)))
setMethod("plot", c(x = "shuffle", y = "missing"),
  function(x, xlab = "u", ylab = "s(u)", ...) {
    m <- length(x@perm)
    bounds <- (0:m) / m

    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
      xlab = xlab, ylab = ylab, ...
    )
    abline(v = bounds, h = bounds, col = "red", lwd = 0.5)

    y0 <- (x@perm - (x@signs > 0)) / m
    y1 <- (x@perm - (x@signs > 0) + x@signs) / m
    segments((0:(m - 1)) / m, y0, (1:m) / m, y1, lwd = 2)
  }
)
