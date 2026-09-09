#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import methods
#' @importFrom graphics abline lines polygon segments
#' @importFrom stats coef cor dbeta integrate pbeta qbeta runif splinefun
## usethis namespace: end
NULL

#' Class of uniform-distribution-preserving transformations
#'
#' `udp` is the virtual superclass of every transformation in the package: a
#' map of the unit interval (or unit hypercube) to itself that sends the
#' uniform distribution to the uniform distribution. \linkS4class{vtransform},
#' \linkS4class{shuffle} and \linkS4class{udpcosine} all inherit from it. The
#' class has no slots and cannot be instantiated; it exists so that
#' [udptrans()] and [udpsi()] can offer a single evaluation and inversion
#' interface for all of them.
#'
#' @seealso [udptrans()] to evaluate one, [udpsi()] to invert one.
#' @export
setClass("udp", contains = "VIRTUAL")

#' Evaluate a uniform-distribution-preserving transformation
#'
#' Applies any \linkS4class{udp} object to values in `[0, 1]`. The construction
#' used depends on the class of `x`: \linkS4class{vtransform},
#' \linkS4class{shuffle} and \linkS4class{udpcosine} objects each carry their
#' own method.
#'
#' @param x an object of class \linkS4class{udp}.
#' @param u a vector, matrix or time series with values in `[0, 1]`.
#'
#' @return An object shaped like `u` with values in `[0, 1]`.
#' @export
#'
#' @examples
#' udptrans(vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
#' udptrans(shuffle(c(2, 1, 3)), c(0, 0.2, 0.5, 0.9, 1))
#' udptrans(udpcosine(2), c(0, 0.25, 0.5, 0.75, 1))
setGeneric("udptrans", function(x, u) standardGeneric("udptrans"))

#' Stochastically invert a uniform-distribution-preserving transformation
#'
#' A [udptrans()] map is in general many-to-one, so it has no ordinary inverse.
#' `udpsi()` returns a single pre-image of `v`, drawn so that the
#' uniform-distribution-preserving property runs backwards too: if `v` is
#' uniform on `[0, 1]` and `Z` is an independent uniform, then `udpsi(x, v, Z)`
#' is uniform on `[0, 1]`. `Z` carries the randomisation.
#'
#' For every \linkS4class{udp} class this is the same computation: [udpinverse()]
#' enumerates the pre-images of each `v` together with their selection
#' probabilities (`1 / |T'|` at each pre-image, normalised over the row), and
#' `udpsi()` picks one pre-image per `v` by inverse-CDF sampling against `Z`.
#' A \linkS4class{shuffle} has a single pre-image and `Z` is ignored;
#' \linkS4class{udpcosine} weights its `degree` pre-images equally;
#' \linkS4class{vtransform} splits between its two branches by the conditional
#' down-probability [vdownprob()].
#'
#' @param x an object of class \linkS4class{udp}.
#' @param v a vector, matrix or time series with values in `[0, 1]`.
#' @param Z a vector of randomisers with values in `[0, 1]`, the same length as
#'   `v`; defaults to a fresh draw from [stats::runif()]. Ignored when `x` is a
#'   \linkS4class{shuffle}.
#' @param ... further arguments forwarded to [udpinverse()]; for
#'   \linkS4class{vtransform} these are the arguments of [vinverse()], such as
#'   `method`, `tol` and `ngrid`.
#'
#' @return An object shaped like `v` with values in `[0, 1]`.
#' @export
#'
#' @examples
#' udpsi(vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
#' udpsi(shuffle(c(2, 1, 3)), c(0.1, 0.5, 0.9))
#' x <- udpcosine(3)
#' set.seed(1)
#' udpsi(x, udptrans(x, runif(5)))
setGeneric("udpsi", function(x, v, Z = runif(length(v)), ...) {
  standardGeneric("udpsi")
})

#' Enumerate the pre-images of a uniform-distribution-preserving transformation
#'
#' Returns every pre-image of each element of `v` under a [udptrans()] map,
#' packed into a matrix. Because the map is in general many-to-one, each `v`
#' may have several pre-images; `udpinverse()` returns them all.
#'
#' The number of pre-images depends on the class of `x`:
#' * \linkS4class{shuffle}: exactly one (the map is a bijection).
#' * \linkS4class{vtransform}: two -- the lower-branch pre-image ([vinverse()])
#'   in column 1 and the upper-branch one in column 2. At `v` equal to `0` or
#'   `1` the branches coincide and both columns hold the same value.
#' * \linkS4class{udpcosine}: up to `degree`, one per linear piece; `v` equal
#'   to `0` or `1` has fewer (the shared troughs and peaks).
#'
#' The result is a numeric matrix with `length(v)` rows and `k` columns, `k`
#' being the largest number of pre-images any `v` can have for this `x`. Row
#' `i` holds the pre-images of `v[i]` in its leading entries, sorted ascending,
#' with `NA` in the trailing entries when `v[i]` has fewer than `k` of them.
#' The count for row `i` is `sum(!is.na(result[i, ]))`.
#'
#' With `prob = TRUE` the result also carries an `n`-by-`k` `"prob"` attribute,
#' aligned column-for-column with the matrix: `attr(., "prob")[i, j]` is the
#' probability with which [udpsi()] selects pre-image `result[i, j]`, namely
#' `1 / |T'|` at that pre-image normalised over the row (`NA` where the
#' pre-image is `NA`). On the measure-zero set where `T'` is undefined at some
#' pre-image the row falls back to equal probabilities.
#'
#' `v` is treated as a plain numeric vector; unlike [udptrans()] and [udpsi()],
#' `udpinverse()` does not copy the attributes of `v` onto its result (the row
#' dimension is `length(v)`, so it cannot).
#'
#' @param x an object of class \linkS4class{udp}.
#' @param v a vector with values in `[0, 1]`.
#' @param prob logical; if `TRUE`, attach the `"prob"` attribute described
#'   above. [udpsi()] calls `udpinverse()` with `prob = TRUE`.
#' @param tol for the \linkS4class{vtransform} method, the convergence tolerance
#'   passed to [vinverse()].
#' @param ... further arguments passed to methods; for \linkS4class{vtransform}
#'   these are the arguments of [vinverse()], such as `method` and `ngrid`.
#'
#' @return A numeric matrix, `length(v)` by `k`, of pre-images: leading entries
#'   filled, trailing entries `NA`, each row sorted ascending; optionally with
#'   a `"prob"` attribute.
#' @export
#'
#' @examples
#' udpinverse(shuffle(c(2, 1, 3)), c(0.2, 0.5, 0.9))
#' udpinverse(vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
#' udpinverse(udpcosine(3), c(0, 0.4, 1))
#' udpinverse(vlinear(0.4), c(0.2, 0.6), prob = TRUE)
setGeneric("udpinverse", function(x, v, prob = FALSE, ...) {
  standardGeneric("udpinverse")
})

#' Probability that stochastic inversion recovers the original value
#'
#' Let `U` be uniform on `[0, 1]` and `V = udptrans(x, U)`. Feeding `V` back
#' through [udpsi()] with an independent uniform randomiser returns *some*
#' pre-image of `V`; `pcoincide()` is the probability that it is `U` itself,
#' `P(udpsi(x, udptrans(x, U)) == U)`.
#'
#' Write `p_1(v), ..., p_k(v)` for the selection probabilities [udpinverse()]
#' attaches to the pre-images of `v`. The uniform-distribution-preserving
#' property makes the true pre-image `U` given `V = v` follow that same
#' distribution, so the chance of drawing it again is the collision probability
#' `sum_j p_j(v)^2`, and
#' \deqn{\mathrm{pcoincide}(x) = \int_0^1 \sum_j p_j(v)^2 \, dv.}
#'
#' The default method evaluates this integral numerically from [udpinverse()].
#' Classes with a closed form override it: a \linkS4class{shuffle} is a
#' bijection so the probability is `1`; a degree-`k` \linkS4class{udpcosine}
#' has `k` equally weighted pre-images so it is `1 / k`; a
#' \linkS4class{vtransform} gives `delta^2 + (1 - delta)^2 + 2 Var(p_down)`.
#'
#' @param x an object of class \linkS4class{udp}.
#'
#' @return A single unnamed probability in `(0, 1]`.
#' @export
#'
#' @examples
#' pcoincide(shuffle(c(2, 1, 3)))
#' pcoincide(udpcosine(4))
#' pcoincide(vlinear(delta = 0.4))
#' pcoincide(udplegendre(3))
setGeneric("pcoincide", function(x) standardGeneric("pcoincide"))

# Numerically integrate sum_j p_j(v)^2 over v in [0, 1], the general formula
# for pcoincide(). `breaks` are interior points of [0, 1] where the integrand
# has a corner (pre-image count changes); integrating piece by piece keeps each
# call to integrate() on a smooth stretch.
integrate_collision <- function(x, breaks = numeric(0)) {
  g <- function(v) {
    P <- attr(udpinverse(x, v, prob = TRUE), "prob")
    rowSums(P^2, na.rm = TRUE)
  }
  cuts <- sort(unique(c(0, breaks[breaks > 0 & breaks < 1], 1)))
  cuts <- cuts[c(TRUE, diff(cuts) > 1e-9)]
  total <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    total <- total + integrate(g, cuts[i], cuts[i + 1L])$value
  }
  total
}

#' @describeIn pcoincide Numerically integrate `sum_j p_j(v)^2` using the
#'   selection probabilities from [udpinverse()]. Serves any \linkS4class{udp}
#'   class without a closed form.
#' @export
setMethod("pcoincide", "udp", function(x) integrate_collision(x))

# Normalise raw per-branch weights into a selection-probability matrix aligned
# with a pre-image matrix. `w` holds 1 / |T'| at each pre-image (any value
# where `present` is FALSE); `present` is `!is.na(<pre-image matrix>)`. Each row
# is scaled to sum to 1. Rows whose weights are non-finite or fail to normalise
# -- the measure-zero set where T' is undefined at some pre-image -- fall back
# to equal probability over that row's pre-images.
finalise_prob <- function(w, present) {
  w[!present] <- 0
  npt <- rowSums(present)
  p <- w / rowSums(w)
  rs <- rowSums(p)
  bad <- !is.finite(rs) | abs(rs - 1) > 1e-6
  if (any(bad)) {
    p[bad, ] <- (present / npt)[bad, ]
  }
  p[!present] <- NA_real_
  dimnames(p) <- NULL
  p
}

#' @describeIn udpsi Draw one pre-image of each `v` by inverse-CDF sampling of
#'   [udpinverse()] against `Z`. This one method serves every
#'   \linkS4class{udp} class.
#' @export
setMethod("udpsi", "udp", function(x, v, Z = runif(length(v)), ...) {
  n <- length(v)
  M <- udpinverse(x, v, prob = TRUE, ...)
  k <- ncol(M)

  if (k == 1L) {
    out <- if (n == 0L) numeric(0) else M[, 1L]
  } else {
    if (length(Z) != n) {
      stop("'Z' must have the same length as 'v'.", call. = FALSE)
    }
    P <- attr(M, "prob")
    P[is.na(P)] <- 0
    cum <- P %*% upper.tri(matrix(0, k, k), diag = TRUE)
    j <- pmin(rowSums(cum < as.numeric(Z)) + 1L, k)
    out <- if (n == 0L) numeric(0) else M[cbind(seq_len(n), j)]
  }

  if (!is.null(attributes(v))) {
    attributes(out) <- attributes(v)
  }
  out
})
