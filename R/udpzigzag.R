# udp transformations built from unequal-width triangle/zigzag waves.
#
# Partition [0, 1] into consecutive pieces of widths w_1, ..., w_n summing to
# 1. If piece i maps affinely and bijectively onto the whole of [0, 1] (slope
# magnitude 1 / w_i), then conditional on U in piece i (probability w_i),
# T(U) is Uniform(0, 1); mixing over i, T(U) is Uniform(0, 1) unconditionally,
# for any choice of positive widths. Continuity forces the pieces to alternate
# direction, giving a generalized (unequal-width) triangle wave. This is the
# complete general piecewise-linear udp family: udpcosine(k) is the equal-
# width case (n = k) and vlinear(delta) is the n = 2 case
# (widths = c(delta, 1 - delta)).
#
# Because every piece maps onto the whole of [0, 1], the raw pre-image weight
# 1 / |T'| at (almost) every pre-image of any v is exactly that piece's width
# -- constant in v, not just constant across pre-images the way udpcosine's
# equal weights are. That collapses pcoincide()'s general integral
# int sum_j p_j(v)^2 dv to the constant sum_j w_j^2.

#' Class of zigzag udp transformations
#'
#' A zigzag udp transformation is the uniform-distribution-preserving map of
#' `[0, 1]` given by an unequal-width triangle wave: `[0, 1]` is partitioned
#' into consecutive pieces, each an affine bijection onto the whole of
#' `[0, 1]`, alternating direction so the map stays continuous. For `U`
#' uniform on `[0, 1]`, `udptrans(x, U)` is again uniform, for any choice of
#' positive piece widths.
#'
#' This is the complete general piecewise-linear udp family:
#' \linkS4class{udpcosine} of degree `k` is the equal-width case, and
#' [vlinear()] is the two-piece case (see [udpzigzag()]).
#'
#' @slot breaks the full, ascending vector of piece boundaries, `c(0, ...,
#'   1)`, length one more than the number of pieces.
#' @slot up logical; does the first piece increase?
#'
#' @seealso [udpzigzag()] to construct one, [udptrans()] to evaluate it.
#' @include udp-package.R
#' @export
#'
#' @references
#' McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
#' of non-monotonic dependence. \href{https://arxiv.org/abs/2512.10828}{arXiv:2512.10828}
setClass("udpzigzag", contains = "udp", slots = list(
  breaks = "numeric", up = "logical"
))

#' Construct a zigzag udp transformation
#'
#' Specify the pieces either by their interior boundaries (`breaks`) or by
#' their widths (`widths`); exactly one of the two must be given.
#'
#' @param breaks the interior breakpoints (`0` and `1` excluded), strictly
#'   increasing and strictly inside `(0, 1)`; an empty vector gives a single
#'   piece. Mutually exclusive with `widths`.
#' @param widths positive piece widths, in order; any positive values are
#'   accepted and normalized to sum to `1`. Mutually exclusive with `breaks`.
#' @param up logical; does the first piece increase? Defaults to `TRUE`.
#'
#' @return An object of class \linkS4class{udpzigzag}.
#' @export
#'
#' @details
#' A single piece (`breaks = numeric(0)` or `widths` of length `1`) reduces
#' to \code{\link{udpid}()} (`up = TRUE`) or \code{\link{udpflip}()}
#' (`up = FALSE`). Equal widths reduce to `udpcosine(n)`, `n` the number of
#' pieces, matching \linkS4class{udpcosine}'s own starting direction with
#' `up = TRUE` for odd `n` and `up = FALSE` for even `n`. Two pieces reduce
#' exactly to `vlinear(delta)` with `widths = c(delta, 1 - delta)`,
#' `up = FALSE`.
#'
#' @examples
#' udpzigzag(breaks = c(0.3, 0.5))
#' udpzigzag(widths = c(1, 2, 1))
#' udpzigzag(widths = c(0.4, 0.6), up = FALSE) # same map as vlinear(0.4)
udpzigzag <- function(breaks = NULL, widths = NULL, up = TRUE) {
  if (is.null(breaks) == is.null(widths)) {
    stop("exactly one of 'breaks' or 'widths' must be given.", call. = FALSE)
  }
  if (!is.logical(up) || length(up) != 1L || is.na(up)) {
    stop("'up' must be a single logical value.", call. = FALSE)
  }
  if (!is.null(widths)) {
    if (!is.numeric(widths) || length(widths) < 1L || anyNA(widths) ||
      any(widths <= 0) || any(!is.finite(widths))) {
      stop("'widths' must be one or more finite positive numbers.",
        call. = FALSE
      )
    }
    full <- c(0, cumsum(widths) / sum(widths))
    full[length(full)] <- 1
  } else {
    if (!is.numeric(breaks) || anyNA(breaks) || any(!is.finite(breaks)) ||
      any(breaks <= 0 | breaks >= 1) ||
      (length(breaks) > 1L && any(diff(breaks) <= 0))) {
      stop(
        "'breaks' must be strictly increasing and strictly inside (0, 1).",
        call. = FALSE
      )
    }
    full <- c(0, breaks, 1)
  }
  new("udpzigzag", breaks = full, up = up)
}

#' @describeIn udptrans Evaluate a zigzag udp transformation.
#' @export
setMethod("udptrans", "udpzigzag", function(x, u) {
  uu <- as.numeric(u)
  b <- x@breaks
  k <- findInterval(uu, b, all.inside = TRUE)
  w <- diff(b)[k]
  increasing <- x@up == ((k %% 2L) == 1L)
  out <- ifelse(increasing, (uu - b[k]) / w, (b[k + 1L] - uu) / w)
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' Roots of a zigzag udp transformation
#'
#' `udptrans()` is an unequal-width triangle wave and so is not injective: a
#' value `v` in `(0, 1)` has one pre-image per piece, while `v = 0` and
#' `v = 1` have fewer (pieces sharing a trough or peak breakpoint collapse to
#' the same root). Each root carries the width of every piece that collapsed
#' onto it, for use as the raw `1 / |T'|` selection weight in
#' [udpinverse()].
#'
#' @param x an object of class \linkS4class{udpzigzag}.
#' @param v a vector with values in `[0, 1]`.
#'
#' @return A list the same length as `v`; element `j` is a list with `u`, the
#'   sorted distinct roots of `udptrans(x, u) == v[j]`, and `w`, the matching
#'   summed piece widths.
#' @keywords internal
udpzigzaginverse <- function(x, v) {
  if (anyNA(v) || any(v < 0 | v > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  b <- x@breaks
  n <- length(b) - 1L
  w <- diff(b)
  k <- seq_len(n)
  increasing <- x@up == ((k %% 2L) == 1L)
  lapply(v, function(vi) {
    u <- ifelse(increasing, b[k] + vi * w, b[k + 1L] - vi * w)
    o <- order(u)
    u <- u[o]
    wt <- w[o]
    grp <- cumsum(c(TRUE, diff(u) > 1e-9))
    list(
      u = unname(vapply(split(u, grp), `[`, 1, FUN.VALUE = 0)),
      w = unname(vapply(split(wt, grp), sum, FUN.VALUE = 0))
    )
  })
}

#' @describeIn udpinverse Pre-images of a zigzag udp transformation: a matrix
#'   with one column per piece holding the roots [udpzigzaginverse()] of each
#'   `v`, sorted ascending and left-packed, `NA`-padded where `v` equal to `0`
#'   or `1` has fewer roots than pieces. With `prob = TRUE` the `"prob"`
#'   attribute weights each root by its piece width (summed over pieces that
#'   collapse onto it), normalized over the row.
#' @export
setMethod("udpinverse", "udpzigzag", function(x, v, prob = FALSE, ...) {
  vv <- as.numeric(v)
  if (anyNA(vv) || any(vv < 0 | vv > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  n <- length(x@breaks) - 1L
  res <- udpzigzaginverse(x, vv)
  roots <- lapply(res, `[[`, "u")
  lens <- lengths(roots)
  M <- matrix(NA_real_, length(vv), n)
  M[cbind(rep(seq_along(roots), lens), sequence(lens))] <-
    unlist(roots, use.names = FALSE)
  if (prob) {
    weights <- lapply(res, `[[`, "w")
    W <- matrix(NA_real_, length(vv), n)
    W[cbind(rep(seq_along(weights), lens), sequence(lens))] <-
      unlist(weights, use.names = FALSE)
    present <- !is.na(M)
    attr(M, "prob") <- finalise_prob(W, present)
  }
  M
})

#' @describeIn udpderiv `1 / w` or `-1 / w`, `w` the width of the piece
#'   containing `u`, the sign matching the direction of that piece, taking
#'   the left piece at a kink and the right piece at `u = 0`, where there is
#'   no left piece.
#' @export
setMethod("udpderiv", "udpzigzag", function(x, u) {
  uu <- as.numeric(u)
  b <- x@breaks
  k <- findInterval(uu - 1e-9, b, all.inside = TRUE)
  w <- diff(b)[k]
  increasing <- x@up == ((k %% 2L) == 1L)
  out <- ifelse(increasing, 1 / w, -1 / w)
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

# Breakpoints: the piece boundaries, where the slope switches sign.
setMethod("udpbreaks", "udpzigzag", function(x) x@breaks)

#' @describeIn pcoincide Every piece maps onto the whole of `[0, 1]`, so each
#'   pre-image's selection probability is its piece's width, constant rather
#'   than varying with `v`; the collision probability `sum_j p_j(v)^2`
#'   reduces to the constant `sum_j w_j^2`, needing no integration.
#' @export
setMethod("pcoincide", "udpzigzag", function(x) sum(diff(x@breaks)^2))

#' Plot method for the udpzigzag class
#'
#' Draws the graph of the zigzag udp transformation by connecting its
#' piece-boundary values directly -- the map is exactly piecewise linear, so
#' no intermediate grid is needed.
#'
#' @param x an object of class \linkS4class{udpzigzag}.
#' @param xlab,ylab axis labels.
#' @param embellish style of the kink gridlines: `"none"` (the default) to
#'   omit them, `"colour"` for red, or `"bw"` for grey.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(udpzigzag(widths = c(1, 2, 1)), embellish = "colour")
#' plot(udpzigzag(breaks = c(0.3, 0.5), up = FALSE), embellish = "bw")
setMethod("plot", c(x = "udpzigzag", y = "missing"),
  function(x, xlab = "u", ylab = "T(u)", embellish = c("none", "colour", "bw"),
           ...) {
    emb <- plot_embellish(embellish)
    b <- x@breaks
    # pty = "s" makes the plotting region square in physical inches; without
    # it, asp = 1 stretches one axis's displayed range past [0, 1] on any
    # device or panel that isn't already exactly square.
    op <- graphics::par(pty = "s")
    on.exit(graphics::par(op))
    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", asp = 1,
      xlab = xlab, ylab = ylab, ...
    )
    if (!is.null(emb)) {
      bu <- b[b > 0 & b < 1]
      if (length(bu)) segments(bu, 0, bu, 1, col = emb$grid, lwd = 0.5)
    }
    lines(b, udptrans(x, b), lwd = 1.5)
  }
)
