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
#' The choice of pre-image depends on the class of `x`:
#' * \linkS4class{vtransform}: the lower pre-image (at or below the fulcrum) is
#'   returned with probability equal to the conditional down-probability
#'   [vdownprob()], the upper one otherwise.
#' * \linkS4class{udpcosine}: one of the `degree` pre-images is returned, root
#'   number `ceiling(Z * degree)` counting from the left.
#' * \linkS4class{shuffle}: a shuffle is a bijection, so its inverse is
#'   deterministic ([shinverse()]) and `Z` is ignored.
#'
#' @param x an object of class \linkS4class{udp}.
#' @param v a vector, matrix or time series with values in `[0, 1]`.
#' @param Z a vector of randomisers with values in `[0, 1]`, the same length as
#'   `v`; defaults to a fresh draw from [stats::runif()]. Ignored when `x` is a
#'   \linkS4class{shuffle}.
#' @param tol for the \linkS4class{vtransform} method, the convergence tolerance
#'   passed to [vinverse()].
#' @param ... further arguments passed to methods; for \linkS4class{vtransform}
#'   these are the arguments of [vinverse()], such as `method`.
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
