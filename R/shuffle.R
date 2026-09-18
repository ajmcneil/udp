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
#' @seealso [shuffle()] to construct one; [udptrans()] and [shinverse()] to
#'   evaluate it and its inverse.
#' @include udp-package.R
#' @export
setClass("shuffle", contains = "udp",
  slots = list(perm = "integer", signs = "numeric"))

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

#' Constructor function for the identity transformation
#'
#' The identity map `T(u) = u`, i.e. `shuffle(1)`: a single strip spanning
#' the whole domain, mapped onto itself with slope `1`.
#'
#' @return An object of class \linkS4class{shuffle}.
#' @export
#'
#' @examples
#' udpid()
#' udptrans(udpid(), c(0, 0.25, 0.5, 0.75, 1))
udpid <- function() {
  shuffle(1)
}

#' Constructor function for the reflection transformation
#'
#' The reflection `T(u) = 1 - u`, i.e. `shuffle(1, signs = -1)`: a single
#' strip spanning the whole domain, mapped onto itself with slope `-1`.
#'
#' @return An object of class \linkS4class{shuffle}.
#' @export
#'
#' @examples
#' udpflip()
#' udptrans(udpflip(), c(0, 0.25, 0.5, 0.75, 1))
udpflip <- function() {
  shuffle(1, signs = -1)
}

#' @describeIn udptrans Evaluate a shuffle.
#' @export
setMethod("udptrans", "shuffle", function(x, u) {
  m <- length(x@perm)
  i <- as.integer(pmax(pmin(floor(u * m) + 1, m), 1))
  t <- u * m - (i - 1)
  s <- x@signs[i]
  out <- (x@perm[i] - (s > 0) + s * t) / m
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' Evaluate the inverse of a shuffle
#'
#' The inverse of a shuffle is itself a shuffle, with permutation
#' `order(perm)` and the signs reordered to match.
#'
#' @param x an object of class \linkS4class{shuffle}.
#' @param v a vector with values in `[0, 1]`.
#'
#' @return An object shaped like `v` with values in `[0, 1]`.
#' @keywords internal
shinverse <- function(x, v) {
  o <- order(x@perm)
  udptrans(new("shuffle", perm = o, signs = x@signs[o]), v)
}

#' @describeIn udpinverse Pre-image of a shuffle: a one-column matrix holding
#'   the single pre-image [shinverse()] of each `v`. With `prob = TRUE` the
#'   `"prob"` attribute is identically `1`.
#' @export
setMethod("udpinverse", "shuffle", function(x, v, prob = FALSE, ...) {
  M <- matrix(shinverse(x, as.numeric(v)), ncol = 1L)
  if (prob) {
    attr(M, "prob") <- finalise_prob(matrix(1, nrow(M), 1L), !is.na(M))
  }
  M
})

#' @describeIn udpderiv The slope `signs[i]` of the strip containing `u`
#'   (`1` or `-1`), taking the left strip at a strip boundary and the right
#'   strip at `u = 0`, where there is no left strip.
#' @export
setMethod("udpderiv", "shuffle", function(x, u) {
  m <- length(x@perm)
  uu <- as.numeric(u)
  i <- as.integer(floor(uu * m - 1e-9)) + 1L
  i <- pmin(pmax(i, 1L), m)
  out <- x@signs[i]
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

# Breakpoints: the m strip boundaries, where the slope can switch sign.
setMethod("udpbreaks", "shuffle", function(x) {
  (0:length(x@perm)) / length(x@perm)
})

#' @describeIn pcoincide A shuffle is a bijection, so stochastic inversion
#'   always recovers the original value: the probability is `1`.
#' @export
setMethod("pcoincide", "shuffle", function(x) 1)

#' Plot method for the shuffle class
#'
#' Draws the graph of the shuffle as one thick black line segment per strip,
#' over thin gridlines at the boundaries of the vertical and horizontal strips.
#'
#' @param x an object of class \linkS4class{shuffle}.
#' @param xlab,ylab axis labels.
#' @param embellish style of the strip-boundary gridlines: `"none"` (the
#'   default) to omit them, `"colour"` for red, or `"bw"` for grey.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(shuffle(c(3, 1, 2), signs = c(1, -1, 1)), embellish = "colour")
#' plot(shuffle(c(3, 1, 2), signs = c(1, -1, 1)), embellish = "none")
setMethod("plot", c(x = "shuffle", y = "missing"),
  function(x, xlab = "u", ylab = "T(u)", embellish = c("none", "colour", "bw"),
           ...) {
    emb <- plot_embellish(embellish)
    m <- length(x@perm)
    bounds <- (0:m) / m

    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", asp = 1,
      xlab = xlab, ylab = ylab, ...
    )
    if (!is.null(emb)) {
      segments(bounds, 0, bounds, 1, col = emb$grid, lwd = 0.5)
      segments(0, bounds, 1, bounds, col = emb$grid, lwd = 0.5)
    }

    y0 <- (x@perm - (x@signs > 0)) / m
    y1 <- (x@perm - (x@signs > 0) + x@signs) / m
    segments((0:(m - 1)) / m, y0, (1:m) / m, y1, lwd = 1.5)
  }
)

# The reflected shuffle u -> 1 - T(u): range strip perm[i] becomes m + 1 -
# perm[i] and every slope flips. udptrans(shuffle_complement(x), u) equals
# 1 - udptrans(x, u) exactly.
shuffle_complement <- function(x) {
  m <- length(x@perm)
  new("shuffle", perm = m + 1L - x@perm, signs = -x@signs)
}

# E[T(U) | U > 1/2] - E[T(U) | U < 1/2] for U uniform on [0, 1]: how much
# higher the shuffle sits on the right half of its domain than on the left. A
# signed measure of upward (> 0) versus downward (< 0) trend that, unlike
# cov(U, T(U)), is not thrown by steep behaviour near the endpoints. T is
# piecewise linear, so its mean over a sub-interval that crosses no strip
# boundary is its value at the sub-interval midpoint; 1/2 is added to the knots
# so no sub-interval straddles it.
shuffle_trend <- function(x) {
  kn <- sort(unique(c((0:length(x@perm)) / length(x@perm), 0.5)))
  mid <- (kn[-1] + kn[-length(kn)]) / 2
  wT <- diff(kn) * udptrans(x, mid)
  2 * (sum(wT[mid > 0.5]) - sum(wT[mid < 0.5]))
}

#' Class of bivariate shuffle fits
#'
#' A `bishuffle` object packages the pair of \linkS4class{shuffle}
#' transformations fitted by [aceshuffle()] to the two columns of a bivariate
#' uniform sample, together with the achieved correlation and the number of
#' alternating sweeps used to find them.
#'
#' @slot shuffle1,shuffle2 the fitted \linkS4class{shuffle} objects for the
#'   two columns.
#' @slot correlation numeric; the achieved linear correlation between
#'   `udptrans(shuffle1, U[, 1])` and `udptrans(shuffle2, U[, 2])`.
#' @slot iterations integer; the number of alternating sweeps performed.
#'
#' @seealso [aceshuffle()] to construct one, [udptrans()] to apply it.
#' @include udp-package.R
#' @export
setClass("bishuffle", slots = list(
  shuffle1 = "shuffle", shuffle2 = "shuffle",
  correlation = "numeric", iterations = "integer"
))

#' @describeIn bishuffle-class Show method for bishuffle objects.
#' @param object an object of class \linkS4class{bishuffle}.
#' @export
setMethod("show", "bishuffle", function(object) {
  cat("An object of class \"bishuffle\"\n")
  cat("m: ", length(object@shuffle1@perm), "\n", sep = "")
  cat("correlation: ", format(object@correlation, digits = 4), "\n", sep = "")
  cat("iterations: ", object@iterations, "\n", sep = "")
})

#' Alternating conditional expectation for shuffle transformations
#'
#' Given a bivariate sample of (approximately) uniform variables,
#' `aceshuffle()` searches for a [shuffle()] of each column that maximizes
#' the linear correlation between the transformed columns.
#'
#' It alternates in the manner of alternating conditional expectations: with the
#' shuffle of column 2 fixed it selects the shuffle of column 1 maximizing the
#' correlation, then with that new shuffle of column 1 fixed it selects the
#' shuffle of column 2, repeating until neither shuffle changes or `maxit`
#' sweeps have been done.
#'
#' Each selection is solved in closed form. Written as a covariance (a shuffle
#' leaves the variance of uniform data essentially unchanged), the objective
#' separates: the sign on each strip is chosen independently, and the
#' permutation follows from the rearrangement inequality. A sweep therefore
#' costs `O(n + m log m)`.
#'
#' The fitted pair is unique only up to the joint reflection
#' `(shuffle1, shuffle2) -> (1 - shuffle1, 1 - shuffle2)`, which leaves the
#' correlation unchanged. `aceshuffle()` returns the orientation whose shuffles
#' trend upward, measured for `U` uniform by
#' `E[shuffle(U) | U > 1/2] - E[shuffle(U) | U < 1/2]` (how much higher the
#' shuffle sits on the right half of its domain than the left). It keeps the
#' orientation with the larger sum of that measure over `shuffle1` and
#' `shuffle2`, breaking a tie (for instance under perfect negative dependence,
#' where the two shuffles cannot both trend upward) on `shuffle1` then
#' `shuffle2`. The measure is a plain trend, not a slope at `u = 1`: it does
#' not force `shuffle1` to match a particular generating transform when that
#' transform is close to symmetric about the centre of the unit square.
#'
#' @param U a two-column numeric matrix with values in `[0, 1]`.
#' @param m the common length of the two permutations (the number of strips).
#' @param maxit maximum number of alternating sweeps.
#' @param init1,init2 optional starting \linkS4class{shuffle} objects with
#'   permutations of length `m`; the identity shuffle is used by default.
#'
#' @return An object of class \linkS4class{bishuffle}.
#' @export
#'
#' @examples
#' set.seed(1)
#' u1 <- runif(2000)
#' u2 <- (u1 + 0.5) %% 1
#' fit <- aceshuffle(cbind(u1, u2), m = 4)
#' fit@correlation
#' plot(fit, embellish = "colour")
aceshuffle <- function(U, m, maxit = 100L, init1 = NULL, init2 = NULL) {
  if (!is.matrix(U) || !is.numeric(U) || ncol(U) != 2L || nrow(U) < 2L) {
    stop("'U' must be a numeric matrix with two columns and at least two rows.",
      call. = FALSE
    )
  }
  U1 <- U[, 1]
  U2 <- U[, 2]
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 1) {
    stop("'m' must be a single positive integer.", call. = FALSE)
  }
  m <- as.integer(m)

  s1 <- if (is.null(init1)) shuffle(seq_len(m)) else init1
  s2 <- if (is.null(init2)) shuffle(seq_len(m)) else init2
  if (!is(s1, "shuffle") || !is(s2, "shuffle") ||
    length(s1@perm) != m || length(s2@perm) != m) {
    stop("'init1' and 'init2' must be shuffles with permutations of length m.",
      call. = FALSE
    )
  }

  # domain strip index (1..m) and within-strip position for each point
  strip_index <- function(u) as.integer(pmax(pmin(floor(u * m) + 1, m), 1))
  i1 <- strip_index(U1)
  i2 <- strip_index(U2)
  t1 <- U1 * m - (i1 - 1)
  t2 <- U2 * m - (i2 - 1)
  f1 <- factor(i1, levels = seq_len(m))
  f2 <- factor(i2, levels = seq_len(m))

  # shuffle of the variable with strip data (fi, t) that best matches target y
  best <- function(fi, t, y) {
    yc <- y - mean(y)
    S <- tapply(yc, fi, sum)
    A <- tapply(t * yc, fi, sum)
    S[is.na(S)] <- 0
    A[is.na(A)] <- 0
    shuffle(rank(S, ties.method = "first"), ifelse(2 * A >= S, 1, -1))
  }

  unchanged <- function(a, b) {
    identical(a@perm, b@perm) && identical(a@signs, b@signs)
  }

  it <- 0L
  repeat {
    it <- it + 1L
    new1 <- best(f1, t1, udptrans(s2, U2))
    new2 <- best(f2, t2, udptrans(new1, U1))
    done <- unchanged(new1, s1) && unchanged(new2, s2)
    s1 <- new1
    s2 <- new2
    if (done || it >= maxit) break
  }

  # resolve the joint reflection (s1, s2) vs (1 - s1, 1 - s2) in favour of the
  # upward-trending orientation; keys compared lexicographically so a zero total
  # trend falls back to s1, then s2
  tr <- c(shuffle_trend(s1), shuffle_trend(s2))
  keys <- c(sum(tr), tr)
  decisive <- keys[abs(keys) > 1e-9]
  if (length(decisive) && decisive[1] < 0) {
    s1 <- shuffle_complement(s1)
    s2 <- shuffle_complement(s2)
  }

  new("bishuffle",
    shuffle1 = s1, shuffle2 = s2,
    correlation = cor(udptrans(s1, U1), udptrans(s2, U2)), iterations = it
  )
}

#' Plot method for the bishuffle class
#'
#' Draws the two shuffles of a fitted \linkS4class{bishuffle} object side by
#' side, in the manner of [graphics::plot()], sharing a common `[0, 1]`
#' frame (every \linkS4class{shuffle} plot already does).
#'
#' @param x an object of class \linkS4class{bishuffle}.
#' @param xlab x-axis label, shared by both panels.
#' @param ylab y-axis label(s) for the `shuffle1` and `shuffle2` panels
#'   respectively; a single value is recycled for both.
#' @param embellish style of the strip-boundary gridlines in each panel,
#'   forwarded to [plot,shuffle,missing-method]: `"none"` (the default),
#'   `"colour"`, or `"bw"`.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' set.seed(1)
#' u1 <- runif(2000)
#' u2 <- (u1 + 0.5) %% 1
#' fit <- aceshuffle(cbind(u1, u2), m = 4)
#' plot(fit, embellish = "colour")
setMethod("plot", c(x = "bishuffle", y = "missing"),
  function(x, xlab = "u", ylab = NULL,
           embellish = c("none", "colour", "bw"), ...) {
    if (is.null(ylab)) {
      ylab <- c("shuffle1(u)", "shuffle2(u)")
    }
    ylab <- rep_len(ylab, 2)
    # pty = "s" makes each side-by-side panel square in physical inches.
    # Without it, the shuffle plot method's asp = 1 forces 1 data-unit to mean
    # the same physical distance on both axes of an oblong panel, which it can
    # only do by stretching one axis's displayed range past [0, 1].
    op <- graphics::par(mfrow = c(1, 2), pty = "s")
    on.exit(graphics::par(op))
    plot(x@shuffle1,
      xlab = xlab, ylab = ylab[1], embellish = embellish,
      main = "shuffle1", ...
    )
    plot(x@shuffle2,
      xlab = xlab, ylab = ylab[2], embellish = embellish,
      main = "shuffle2", ...
    )
  }
)

#' @describeIn udptrans Apply the pair of shuffles in a \linkS4class{bishuffle}
#'   object to the two columns of a bivariate matrix, `cbind(udptrans(shuffle1,
#'   u[, 1]), udptrans(shuffle2, u[, 2]))`. Unlike a single
#'   \linkS4class{shuffle}, which acts on a vector, `u` here must be a
#'   two-column matrix since the two columns get different shuffles.
#' @export
setMethod("udptrans", "bishuffle", function(x, u) {
  if (!is.matrix(u) || ncol(u) != 2L) {
    stop("'u' must be a two-column matrix.", call. = FALSE)
  }
  cbind(udptrans(x@shuffle1, u[, 1]), udptrans(x@shuffle2, u[, 2]))
})
