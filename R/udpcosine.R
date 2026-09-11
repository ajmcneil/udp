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
#'
#' @references
#' McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
#' of non-monotonic dependence. \href{https://arxiv.org/abs/2512.10828}{arXiv:2512.10828}
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
#' @keywords internal
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

#' @describeIn udpinverse Pre-images of a cosine udp transformation: a matrix
#'   with `degree` columns holding the roots [udpcosinverse()] of each `v`,
#'   sorted ascending and left-packed, `NA`-padded where `v` equal to `0` or
#'   `1` has fewer than `degree` roots. With `prob = TRUE` the `"prob"`
#'   attribute is equal over each row's roots.
#' @export
setMethod("udpinverse", "udpcosine", function(x, v, prob = FALSE, ...) {
  vv <- as.numeric(v)
  if (anyNA(vv) || any(vv < 0 | vv > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  k <- x@degree
  roots <- udpcosinverse(x, vv)
  lens <- lengths(roots)
  M <- matrix(NA_real_, length(vv), k)
  M[cbind(rep(seq_along(roots), lens), sequence(lens))] <-
    unlist(roots, use.names = FALSE)
  if (prob) {
    present <- !is.na(M)
    attr(M, "prob") <- finalise_prob(present + 0, present)
  }
  M
})

#' @describeIn udpderiv `degree` or `-degree`, the sign matching the direction
#'   of the linear piece containing `u`, taking the left piece at a kink and
#'   the right piece at `u = 0`, where there is no left piece.
#' @export
setMethod("udpderiv", "udpcosine", function(x, u) {
  d <- x@degree
  uu <- as.numeric(u)
  k <- as.integer(floor(uu * d - 1e-9)) + 1L
  k <- pmin(pmax(k, 1L), d)
  increasing <- (k %% 2L) == (d %% 2L)
  out <- ifelse(increasing, d, -d)
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' @describeIn pcoincide The degree-`k` triangle wave has `k` equally weighted
#'   pre-images at almost every `v`, so the probability is `1 / k`.
#' @export
setMethod("pcoincide", "udpcosine", function(x) 1 / x@degree)

#' Plot method for the udpcosine class
#'
#' Draws the graph of the cosine udp transformation over thin gridlines at the
#' boundaries of its `degree` linear pieces.
#'
#' @param x an object of class \linkS4class{udpcosine}.
#' @param xlab,ylab axis labels.
#' @param embellish style of the kink gridlines: `"colour"` (the default) for
#'   red, `"bw"` for grey, or `"none"` to omit them.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(udpcosine(3))
#' plot(udpcosine(4), embellish = "bw")
setMethod("plot", c(x = "udpcosine", y = "missing"),
  function(x, xlab = "u", ylab = "T(u)", embellish = c("colour", "bw", "none"),
           ...) {
    emb <- plot_embellish(embellish)
    u <- (0:x@degree) / x@degree
    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
      xlab = xlab, ylab = ylab, ...
    )
    if (!is.null(emb)) segments(u, 0, u, 1, col = emb$grid, lwd = 0.5)
    lines(u, udptrans(x, u), lwd = 1.5)
  }
)
