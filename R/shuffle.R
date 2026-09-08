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

#' Alternating conditional expectation for shuffle transformations
#'
#' Given a paired sample of (approximately) uniform variables, `aceshuffle()`
#' searches for a [shuffle()] of each margin that maximises the linear
#' correlation between the transformed variables.
#'
#' It alternates in the manner of alternating conditional expectations: with the
#' shuffle of `U2` fixed it selects the shuffle of `U1` maximising the
#' correlation, then with that new shuffle of `U1` fixed it selects the shuffle
#' of `U2`, repeating until neither shuffle changes or `maxit` sweeps have been
#' done.
#'
#' Each selection is solved in closed form. Written as a covariance (a shuffle
#' leaves the variance of uniform data essentially unchanged), the objective
#' separates: the sign on each strip is chosen independently, and the
#' permutation follows from the rearrangement inequality. A sweep therefore
#' costs `O(n + m log m)`.
#'
#' @param U1,U2 numeric vectors of equal length with values in `[0, 1]`.
#' @param m the common length of the two permutations (the number of strips).
#' @param maxit maximum number of alternating sweeps.
#' @param init1,init2 optional starting \linkS4class{shuffle} objects with
#'   permutations of length `m`; the identity shuffle is used by default.
#'
#' @return A list with elements `data` (the `cbind(U1, U2)` matrix),
#'   `shuffle1` and `shuffle2` (the fitted \linkS4class{shuffle} objects),
#'   `correlation` (the achieved linear correlation) and `iterations` (the
#'   number of sweeps performed).
#' @export
#'
#' @examples
#' set.seed(1)
#' u1 <- runif(2000)
#' u2 <- (u1 + 0.5) %% 1
#' fit <- aceshuffle(u1, u2, m = 4)
#' fit$correlation
#' plot(fit$shuffle1)
aceshuffle <- function(U1, U2, m, maxit = 100L, init1 = NULL, init2 = NULL) {
  n <- length(U1)
  if (!is.numeric(U1) || !is.numeric(U2) || length(U2) != n || n < 2L) {
    stop("'U1' and 'U2' must be numeric vectors of the same length (>= 2).",
      call. = FALSE
    )
  }
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
    new1 <- best(f1, t1, shtrans(s2, U2))
    new2 <- best(f2, t2, shtrans(new1, U1))
    done <- unchanged(new1, s1) && unchanged(new2, s2)
    s1 <- new1
    s2 <- new2
    if (done || it >= maxit) break
  }

  list(
    data = cbind(U1 = U1, U2 = U2),
    shuffle1 = s1,
    shuffle2 = s2,
    correlation = cor(shtrans(s1, U1), shtrans(s2, U2)),
    iterations = it
  )
}
