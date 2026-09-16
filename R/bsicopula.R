# Bivariate stochastic inversion copulas: couple two udp transformations
# through a copula on their "carrier" uniforms (V1, V2), optionally with
# extra dependence injected via the randomizers udpsi() uses to pick among
# each transformation's pre-images. 'copula' and 'rvinecopulib' are Suggests,
# not Imports, so every reference to their classes/functions below goes
# through requireNamespace()/inherits() and the :: operator rather than a
# formal S4 class union or NAMESPACE import.

# TRUE if x is a parCopula object from the copula package -- FALSE, not an
# error, when copula isn't installed.
is_parCopula <- function(x) {
  requireNamespace("copula", quietly = TRUE) && methods::is(x, "parCopula")
}

# TRUE if x is a bicop_dist object from rvinecopulib. inherits() only looks
# at the class attribute, so this needs no namespace load.
is_bicop_dist <- function(x) inherits(x, "bicop_dist")

# Sample n draws from a base copula that is either a parCopula or a
# bicop_dist object, dispatching to the matching package's own sampler.
sample_basecopula <- function(n, basecopula) {
  if (is_parCopula(basecopula)) {
    copula::rCopula(n, basecopula)
  } else {
    rvinecopulib::rbicop(n, basecopula)
  }
}

#' Class of simplified D-vine randomizer models
#'
#' A randomizer model for \linkS4class{bsicopula}: the three non-trivial
#' pair-copulas of a 4-dimensional D-vine on `(Z1, V1, V2, Z2)`, whose two
#' outer tree-1 edges, `C_{Z1,V1}` and `C_{V2,Z2}`, are fixed to the
#' independence copula and are not slots of this class -- exactly the
#' constraint stochastic inversion requires, since [udpsi()]'s randomizer
#' must stay independent of the value it randomizes. `copZ1Z2_V1V2` is the
#' tree-3 edge; `copZ1V2_V1` and `copV1Z2_V2` are the two remaining tree-2
#' edges.
#'
#' @slot copZ1Z2_V1V2 a bicop_dist object (\pkg{rvinecopulib}), the tree-3
#'   pair-copula of `Z1` and `Z2` given `(V1, V2)`.
#' @slot copZ1V2_V1 a bicop_dist object, the tree-2 pair-copula of `Z1` and
#'   `V2` given `V1`.
#' @slot copV1Z2_V2 a bicop_dist object, the tree-2 pair-copula of `V1` and
#'   `Z2` given `V2`.
#'
#' @seealso [sdvine()] to construct one; \linkS4class{bsicopula} to use it.
#' @export
setClass("sdvine", slots = list(
  copZ1Z2_V1V2 = "ANY",
  copZ1V2_V1   = "ANY",
  copV1Z2_V2   = "ANY"
))

#' Construct a simplified D-vine randomizer model
#'
#' @param copZ1Z2_V1V2 a bicop_dist object (\pkg{rvinecopulib}), the tree-3
#'   pair-copula of `Z1` and `Z2` given `(V1, V2)`. Must be supplied.
#' @param copZ1V2_V1,copV1Z2_V2 bicop_dist objects, the two tree-2
#'   pair-copulas. Default to `rvinecopulib::bicop_dist()`, the independence
#'   copula.
#'
#' @return An object of class \linkS4class{sdvine}.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   sdvine(rvinecopulib::bicop_dist("gaussian", 0, 0.4))
#' }
sdvine <- function(copZ1Z2_V1V2,
                    copZ1V2_V1 = rvinecopulib::bicop_dist(),
                    copV1Z2_V2 = rvinecopulib::bicop_dist()) {
  if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("Package 'rvinecopulib' is required to construct an 'sdvine' object.",
      call. = FALSE
    )
  }
  if (!is_bicop_dist(copZ1Z2_V1V2) || !is_bicop_dist(copZ1V2_V1) || !is_bicop_dist(copV1Z2_V2)) {
    stop(
      "'copZ1Z2_V1V2', 'copZ1V2_V1' and 'copV1Z2_V2' must all be bicop_dist objects.",
      call. = FALSE
    )
  }
  new("sdvine",
    copZ1Z2_V1V2 = copZ1Z2_V1V2, copZ1V2_V1 = copZ1V2_V1, copV1Z2_V2 = copV1Z2_V2
  )
}

#' Class of bivariate stochastic inversion copulas
#'
#' Couples two \linkS4class{udp} transformations through a copula on their
#' "carrier" uniforms `(V1, V2)`. Margin `i`'s draw is `udpsi(udp_i, V_i,
#' Z_i)`: [udpsi()] applied to the base draw `V_i`, with the pre-image
#' selected by the randomizer `Z_i`. `Z_i` defaults to an independent
#' uniform, giving margins that agree with `basecopula` up to the
#' pre-image-selection randomness. When `randomizermod` is instead an
#' \linkS4class{sdvine} model, `(Z1, Z2)` are drawn dependently -- on `(V1,
#' V2)` and on each other -- injecting further dependence between the two
#' margins, including non-monotonic dependence when `udp1`/`udp2` are
#' many-to-one.
#'
#' @slot basecopula a parCopula object (\pkg{copula}) or a bicop_dist object
#'   (\pkg{rvinecopulib}), the copula of `(V1, V2)`.
#' @slot udp1,udp2 objects of class \linkS4class{udp}, applied via
#'   stochastic inversion ([udpsi()]) to `V1` and `V2` respectively.
#' @slot randomizermod `NULL` (independent randomizers) or an
#'   \linkS4class{sdvine} object.
#'
#' @seealso [bsicopula()] to construct one; [rbsicopula()] to sample from one.
#' @include udp-package.R
#' @export
setClass("bsicopula", slots = list(
  basecopula = "ANY",
  udp1 = "udp",
  udp2 = "udp",
  randomizermod = "ANY"
))

#' Construct a bivariate stochastic inversion copula
#'
#' @param basecopula a parCopula object (\pkg{copula}) or a bicop_dist
#'   object (\pkg{rvinecopulib}), the copula of the two transformations'
#'   carrier uniforms `(V1, V2)`.
#' @param udp1,udp2 objects of class \linkS4class{udp}.
#' @param randomizermod `NULL` (the default; independent randomizers) or an
#'   \linkS4class{sdvine} object. When given, `basecopula` must be a
#'   bicop_dist object -- the D-vine sampling machinery uses
#'   \pkg{rvinecopulib}'s h-functions throughout.
#'
#' @return An object of class \linkS4class{bsicopula}.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), udpcosine(2), udpcosine(3))
#' }
bsicopula <- function(basecopula, udp1, udp2, randomizermod = NULL) {
  if (!is_parCopula(basecopula) && !is_bicop_dist(basecopula)) {
    stop(
      "'basecopula' must be a parCopula object (copula package) or a bicop_dist object (rvinecopulib).",
      call. = FALSE
    )
  }
  if (!methods::is(udp1, "udp") || !methods::is(udp2, "udp")) {
    stop("'udp1' and 'udp2' must be objects of class 'udp'.", call. = FALSE)
  }
  if (!is.null(randomizermod)) {
    if (!methods::is(randomizermod, "sdvine")) {
      stop("'randomizermod' must be NULL or an object of class 'sdvine'.", call. = FALSE)
    }
    if (!is_bicop_dist(basecopula)) {
      stop(
        "'basecopula' must be a bicop_dist object when 'randomizermod' is an 'sdvine'.",
        call. = FALSE
      )
    }
  }
  new("bsicopula",
    basecopula = basecopula, udp1 = udp1, udp2 = udp2, randomizermod = randomizermod
  )
}

# (Z1, Z2) given (V1, V2) from the simplified D-vine: basecopula is
# C_{V1,V2} (tree 1), randomizermod's three slots are the tree-2/tree-3
# edges; the two tree-1 outer edges C_{Z1,V1} and C_{V2,Z2} are the
# independence copula, not stored anywhere (see sdvine()). Uses the
# rvinecopulib h-function convention throughout: hbicop(cbind(cond, target),
# cond_var = 1, family = bicop) is C(target | cond), and inverse = TRUE
# solves it for target given a probability level.
sdvine_sample <- function(V1, V2, basecopula, randomizermod) {
  n <- length(V1)
  w1 <- runif(n)
  w2 <- runif(n)

  e21 <- rvinecopulib::hbicop(cbind(V1, V2), cond_var = 1, family = basecopula)
  e12 <- rvinecopulib::hbicop(cbind(V2, V1), cond_var = 1, family = basecopula)

  Z1 <- rvinecopulib::hbicop(cbind(e21, w1),
    cond_var = 1, family = randomizermod@copZ1V2_V1, inverse = TRUE
  )
  e1 <- w1

  y <- rvinecopulib::hbicop(cbind(e1, w2),
    cond_var = 1, family = randomizermod@copZ1Z2_V1V2, inverse = TRUE
  )
  Z2 <- rvinecopulib::hbicop(cbind(e12, y),
    cond_var = 1, family = randomizermod@copV1Z2_V2, inverse = TRUE
  )

  cbind(Z1 = Z1, Z2 = Z2)
}

#' Random sample from a bivariate stochastic inversion copula
#'
#' 1. Draws `(V1, V2)` from `basecopula`.
#' 2. If `randomizermod` is `NULL`, returns `(udpsi(udp1, V1), udpsi(udp2,
#'    V2))` -- each with its own fresh, independent randomizer -- and stops.
#' 3. Otherwise draws `(Z1, Z2)` given `(V1, V2)` from the simplified D-vine
#'    specified by `basecopula` (as `C_{V1,V2}`) and `randomizermod`.
#' 4. Returns `(udpsi(udp1, V1, Z1), udpsi(udp2, V2, Z2))`.
#'
#' @param n number of draws.
#' @param object an object of class \linkS4class{bsicopula}.
#'
#' @return A numeric matrix with `n` rows and columns `U1`, `U2`, each
#'   marginally uniform on `[0, 1]`.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   set.seed(1)
#'   bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), udpcosine(2), udpcosine(3))
#'   rbsicopula(1000, bc)
#' }
rbsicopula <- function(n, object) {
  if (!methods::is(object, "bsicopula")) {
    stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
  }
  V <- sample_basecopula(n, object@basecopula)
  V1 <- V[, 1]
  V2 <- V[, 2]

  if (is.null(object@randomizermod)) {
    U1 <- udpsi(object@udp1, V1)
    U2 <- udpsi(object@udp2, V2)
  } else {
    Z <- sdvine_sample(V1, V2, object@basecopula, object@randomizermod)
    U1 <- udpsi(object@udp1, V1, Z[, "Z1"])
    U2 <- udpsi(object@udp2, V2, Z[, "Z2"])
  }
  cbind(U1 = U1, U2 = U2)
}
