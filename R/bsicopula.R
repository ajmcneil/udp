# Bivariate stochastic inversion copulas: couple two udp transformations
# through a copula on their "carrier" uniforms (V1, V2), optionally with
# extra dependence injected via the randomizers udpsi() uses to pick among
# each transformation's pre-images. 'copula' and 'rvinecopulib' are Suggests,
# not Imports, so every reference to their classes/functions below goes
# through requireNamespace()/inherits() and the :: operator rather than a
# formal S4 class union or NAMESPACE import.

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
#' Each pair-copula is specified in the order its variables appear in the
#' D-vine `(Z1, V1, V2, Z2)`: the first argument of `copZ1Z2_V1V2` is `Z1`
#' and the second `Z2`; of `copZ1V2_V1`, `Z1` then `V2`; of `copV1Z2_V2`,
#' `V1` then `Z2`. This matters only for copulas that are not exchangeable,
#' such as the 90 and 270 degree rotations: transposing one swaps the two.
#'
#' @slot copZ1Z2_V1V2 a bicop_dist object (\pkg{rvinecopulib}), the tree-3
#'   pair-copula of `(Z1, Z2)` given `(V1, V2)`, in that order.
#' @slot copZ1V2_V1 a bicop_dist object, the tree-2 pair-copula of `(Z1, V2)`
#'   given `V1`, in that order.
#' @slot copV1Z2_V2 a bicop_dist object, the tree-2 pair-copula of `(V1, Z2)`
#'   given `V2`, in that order.
#'
#' @seealso [randsdvine()] to construct one; \linkS4class{bsicopula} to use it.
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
#' @export
setClass("randsdvine", slots = list(
  copZ1Z2_V1V2 = "ANY",
  copZ1V2_V1   = "ANY",
  copV1Z2_V2   = "ANY"
))

#' Construct a simplified D-vine randomizer model
#'
#' @param copZ1Z2_V1V2 a bicop_dist object (\pkg{rvinecopulib}), the tree-3
#'   pair-copula of `(Z1, Z2)` given `(V1, V2)`, in that order. Must be
#'   supplied.
#' @param copZ1V2_V1,copV1Z2_V2 bicop_dist objects, the two tree-2
#'   pair-copulas: of `(Z1, V2)` given `V1`, and of `(V1, Z2)` given `V2`, each
#'   in the order shown, which is the order of the variables in the D-vine
#'   `(Z1, V1, V2, Z2)`. The order matters for copulas that are not
#'   exchangeable, such as the 90 and 270 degree rotations. Default to
#'   `rvinecopulib::bicop_dist()`, the independence copula.
#'
#' @return An object of class \linkS4class{randsdvine}.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   randsdvine(rvinecopulib::bicop_dist("gaussian", 0, 0.4))
#' }
randsdvine <- function(copZ1Z2_V1V2,
                    copZ1V2_V1 = rvinecopulib::bicop_dist(),
                    copV1Z2_V2 = rvinecopulib::bicop_dist()) {
  if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("Package 'rvinecopulib' is required to construct a 'randsdvine' object.",
      call. = FALSE
    )
  }
  if (!is_bicop_dist(copZ1Z2_V1V2) || !is_bicop_dist(copZ1V2_V1) || !is_bicop_dist(copV1Z2_V2)) {
    stop(
      "'copZ1Z2_V1V2', 'copZ1V2_V1' and 'copV1Z2_V2' must all be bicop_dist objects.",
      call. = FALSE
    )
  }
  new("randsdvine",
    copZ1Z2_V1V2 = copZ1Z2_V1V2, copZ1V2_V1 = copZ1V2_V1, copV1Z2_V2 = copV1Z2_V2
  )
}

#' Class of mixture-of-copulas randomizer models
#'
#' A randomizer model for \linkS4class{bsicopula} that draws `(Z1, Z2)`
#' given `(V1, V2)` from one of two copulas, `cop1` or `cop2`, chosen by a
#' user-supplied `selector(v1, v2)`. Unlike \linkS4class{randsdvine}, this needs
#' no constraint on `selector` to keep `udpsi()`'s stochastic-inversion
#' guarantee valid: whichever of `cop1`/`cop2` gets picked, it is still some
#' copula, and every copula's own coordinate margins are uniform on
#' `[0, 1]` by definition regardless of the other coordinate or the
#' dependence parameter -- so `Z1 | V1 = v1, V2 = v2` is uniform on
#' `[0, 1]` for every `(v1, v2)`, whatever `selector` does, which is
#' exactly what `Z1` independent of `V1` requires.
#'
#' `selector` must nonetheless be a pure, deterministic function of
#' `(v1, v2)`: [rbsicopula()] calls it once to sample, while
#' [dbsicopula()]/`plot()` call it again, independently, to evaluate the
#' density, and [dbsicopula()] itself calls it once for every corner of an
#' inclusion-exclusion rectangle. If `selector` has any internal randomness
#' or depends on mutable external state, these calls can disagree, silently
#' corrupting results -- so avoid things like
#' `function(v1, v2) pmax(v1, v2) > runif(1)`, and capture any threshold in
#' the closure by value, not by a variable that might be reassigned later.
#' `selector` must also be vectorized (`pmax()`/`ifelse()`/`&`/`|`, not
#' `if`/`&&`/`||`) and must never return `NA`.
#'
#' @slot cop1,cop2 parCopula objects (\pkg{copula}), bicop_dist objects
#'   (\pkg{rvinecopulib}) or \linkS4class{astcopula} objects, the two copulas
#'   of `(Z1, Z2)` to choose between.
#'   `cop1` and `cop2` need not share a backend -- one may be a parCopula
#'   and the other a bicop_dist object -- since sampling and density/CDF
#'   evaluation dispatch on each of `cop1`/`cop2` individually.
#' @slot selector a function `selector(v1, v2)` returning a logical vector
#'   the same length as `v1`/`v2`, with no `NA`s: `TRUE` selects `cop1`,
#'   `FALSE` selects `cop2`. Any further parameters (such as a threshold)
#'   should be captured in `selector`'s closure rather than passed
#'   separately. See Details for the purity requirement this implies.
#'
#' @seealso [randmixture()] to construct one; \linkS4class{bsicopula} to use it.
#' @export
setClass("randmixture", slots = list(
  cop1 = "ANY",
  cop2 = "ANY",
  selector = "function"
))

#' Construct a mixture-of-copulas randomizer model
#'
#' @param cop1,cop2 parCopula objects (\pkg{copula}), bicop_dist objects
#'   (\pkg{rvinecopulib}) or \linkS4class{astcopula} objects, the two copulas
#'   of `(Z1, Z2)` to choose between.
#' @param selector a function `selector(v1, v2)` returning a logical vector
#'   the same length as `v1`/`v2`, with no `NA`s: `TRUE` selects `cop1`,
#'   `FALSE` selects `cop2`. Capture any further parameters in its closure,
#'   e.g. `function(v1, v2) pmax(v1, v2) > 0.7`. Must be pure and
#'   vectorized -- see \linkS4class{randmixture}'s Details.
#'
#' @return An object of class \linkS4class{randmixture}.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   randmixture(
#'     rvinecopulib::bicop_dist("gaussian", 0, 1),
#'     rvinecopulib::bicop_dist("gaussian", 0, -1),
#'     selector = function(v1, v2) pmax(v1, v2) > 0.7
#'   )
#' }
randmixture <- function(cop1, cop2, selector) {
  if (!basecopula_supported(cop1) || !basecopula_supported(cop2)) {
    stop(
      "'cop1' and 'cop2' must each be a parCopula object (copula package), a bicop_dist object (rvinecopulib) or another supported base copula.",
      call. = FALSE
    )
  }
  if (!is.function(selector)) {
    stop("'selector' must be a function.", call. = FALSE)
  }
  new("randmixture", cop1 = cop1, cop2 = cop2, selector = selector)
}

#' Class of bivariate stochastic inversion copulas
#'
#' Couples two \linkS4class{udp} transformations through a copula on their
#' "carrier" uniforms `(V1, V2)`. Margin `i`'s draw is `udpsi(udp_i, V_i,
#' Z_i)`: [udpsi()] applied to the base draw `V_i`, with the pre-image
#' selected by the randomizer `Z_i`. `Z_i` defaults to an independent
#' uniform, giving margins that agree with `basecopula` up to the
#' pre-image-selection randomness. When `randomizermod` is instead an
#' \linkS4class{randsdvine} or \linkS4class{randmixture} model, `(Z1, Z2)` are
#' drawn dependently -- on `(V1, V2)` and on each other -- injecting
#' further dependence between the two margins, including non-monotonic
#' dependence when `udp1`/`udp2` are many-to-one.
#'
#' **Base copulas.** The copula of `(V1, V2)` can be a parametric family from
#' \pkg{rvinecopulib} (`bicop_dist`), a `parCopula` from \pkg{copula}, or the
#' absolute spherical t copula, \linkS4class{astcopula}, the copula of the
#' absolute values of a bivariate t vector with correlation zero (Dias, Han
#' and McNeil, 2027), which has upper tail dependence and is the copula the
#' symmetric v-transform produces from the t copula. See [astcopula()] for its
#' density, CDF, h-functions, Kendall's tau and the calibration of its degrees
#' of freedom to a Kendall's tau. Which base copulas each function supports is
#' stated in its help page; in particular [fitbsicopula()] estimates
#' `bicop_dist` and `astcopula` base copulas.
#'
#' @slot basecopula a parCopula object (\pkg{copula}), a bicop_dist object
#'   (\pkg{rvinecopulib}) or an \linkS4class{astcopula}, the copula of `(V1, V2)`.
#' @slot udp1,udp2 objects of class \linkS4class{udp}, applied via
#'   stochastic inversion ([udpsi()]) to `V1` and `V2` respectively.
#' @slot randomizermod `NULL` (independent randomizers), an
#'   \linkS4class{randsdvine} object, or a \linkS4class{randmixture} object.
#'
#' @seealso [bsicopula()] to construct one; [rbsicopula()] to sample from one.
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
#'
#' Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#' D-vines for stochastic volatility. *Journal of Multivariate Analysis*,
#' **217**, 105695. \doi{10.1016/j.jmva.2026.105695}
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
#' @param basecopula a parCopula object (\pkg{copula}), a bicop_dist
#'   object (\pkg{rvinecopulib}) or an \linkS4class{astcopula}, the copula of the two transformations'
#'   carrier uniforms `(V1, V2)`.
#' @param udp1,udp2 objects of class \linkS4class{udp}.
#' @param randomizermod `NULL` (the default; independent randomizers), an
#'   \linkS4class{randsdvine} object, or a \linkS4class{randmixture} object. A
#'   \linkS4class{randsdvine} needs the h-functions of the base copula, which
#'   any base copula supplies; for a parCopula they are central differences of
#'   its CDF (accurate to about `1e-9` but slower), for a bicop_dist and an
#'   \linkS4class{astcopula} they are exact. The copulas inside a
#'   \linkS4class{randsdvine} must be bicop_dist objects.
#'
#' @return An object of class \linkS4class{bsicopula}.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), udpcosine(2), udpcosine(3))
#' }
bsicopula <- function(basecopula, udp1, udp2, randomizermod = NULL) {
  if (!basecopula_supported(basecopula)) {
    stop(
      "'basecopula' must be a parCopula object (copula package), a bicop_dist object (rvinecopulib) or another supported base copula.",
      call. = FALSE
    )
  }
  if (!methods::is(udp1, "udp") || !methods::is(udp2, "udp")) {
    stop("'udp1' and 'udp2' must be objects of class 'udp'.", call. = FALSE)
  }
  if (!is.null(randomizermod)) {
    if (!methods::is(randomizermod, "randsdvine") && !methods::is(randomizermod, "randmixture")) {
      stop(
        "'randomizermod' must be NULL, an object of class 'randsdvine', or an object of class 'randmixture'.",
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
# independence copula, not stored anywhere (see randsdvine()).
#
# rvinecopulib h-function convention: for (A, B) ~ bicop, hbicop(cbind(a, b),
# cond_var = 1) is P(B <= b | A = a) and cond_var = 2 is P(A <= a | B = b);
# with inverse = TRUE, cond_var = 1 takes cbind(a, level) and returns b, and
# cond_var = 2 takes cbind(level, b) and returns a. Every copula is used in
# the order its variables appear in the D-vine (Z1, V1, V2, Z2) -- base
# (V1, V2), copZ1V2_V1 (Z1, V2), copV1Z2_V2 (V1, Z2), copZ1Z2_V1V2
# (Z1, Z2) -- and the arguments are never swapped: transposing a copula
# that is not exchangeable (a 90 or 270 degree rotation) gives another
# copula. So:
#   e21 = P(V2 <= v2 | V1 = v1)            basecopula_cond_cdf(base, v1, v2, 1)
#   e12 = P(V1 <= v1 | V2 = v2)            basecopula_cond_cdf(base, v2, v1, 2)
#   Z1 | (V1, V2): P(Z1 <= z1 | V2-level e21) under copZ1V2_V1 (Z1, V2),
#                  conditioning on the second variable: cond_var = 2
#   Z2 | (V1, V2, Z1): tree 3 conditions on the first variable, Z1
#   Z2 from P(Z2 <= z2 | V1-level e12) under copV1Z2_V2 (V1, Z2),
#                  conditioning on the first variable: cond_var = 1
randsdvine_sample <- function(V1, V2, basecopula, randomizermod) {
  n <- length(V1)
  w1 <- runif(n)
  w2 <- runif(n)

  e21 <- basecopula_cond_cdf(basecopula, V1, V2, 1L)
  e12 <- basecopula_cond_cdf(basecopula, V2, V1, 2L)

  Z1 <- rvinecopulib::hbicop(cbind(w1, e21),
    cond_var = 2, family = randomizermod@copZ1V2_V1, inverse = TRUE
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

# Checks a randmixture selector(v1, v2) result: must be a logical vector,
# the same length as v1/v2, with no NA (an NA would otherwise reach
# Z1[sel] <- .../Z1[!sel] <- ... in randmixture_sample() and fail there with
# a cryptic "NAs are not allowed in subscripted assignments" instead of
# pointing at the real cause).
validate_selector_result <- function(sel, n) {
  if (!is.logical(sel) || length(sel) != n) {
    stop(
      "'selector' must return a logical vector the same length as 'v1'/'v2'.",
      call. = FALSE
    )
  }
  if (anyNA(sel)) {
    stop("'selector' must not return NA.", call. = FALSE)
  }
}

# (Z1, Z2) given (V1, V2) from a randmixture model: evaluate selector(V1, V2)
# on the realized pair, then draw from cop1 where TRUE and cop2 where FALSE.
# Z1 independent of V1 (and Z2 of V2) holds for any selector -- see
# randmixture-class's documentation -- so unlike randsdvine_sample() this needs no
# h-function machinery, just an ordinary draw from whichever copula was
# picked.
randmixture_sample <- function(V1, V2, randomizermod) {
  n <- length(V1)
  sel <- randomizermod@selector(V1, V2)
  validate_selector_result(sel, n)
  Z1 <- numeric(n)
  Z2 <- numeric(n)
  if (any(sel)) {
    Z <- basecopula_sample(randomizermod@cop1, sum(sel))
    Z1[sel] <- Z[, 1]
    Z2[sel] <- Z[, 2]
  }
  if (any(!sel)) {
    Z <- basecopula_sample(randomizermod@cop2, sum(!sel))
    Z1[!sel] <- Z[, 1]
    Z2[!sel] <- Z[, 2]
  }
  cbind(Z1 = Z1, Z2 = Z2)
}

#' Random sample from a bivariate stochastic inversion copula
#'
#' 1. Draws `(V1, V2)` from `basecopula`.
#' 2. If `randomizermod` is `NULL`, returns `(udpsi(udp1, V1), udpsi(udp2,
#'    V2))` -- each with its own fresh, independent randomizer -- and stops.
#' 3. Otherwise draws `(Z1, Z2)` given `(V1, V2)`: from the simplified
#'    D-vine specified by `basecopula` (as `C_{V1,V2}`) and `randomizermod`
#'    when it is a \linkS4class{randsdvine}, or from `randomizermod@cop1` /
#'    `cop2` according to `randomizermod@selector(V1, V2)` when it is a
#'    \linkS4class{randmixture}.
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
  V <- basecopula_sample(object@basecopula, n)
  V1 <- V[, 1]
  V2 <- V[, 2]

  if (is.null(object@randomizermod)) {
    U1 <- udpsi(object@udp1, V1)
    U2 <- udpsi(object@udp2, V2)
  } else {
    Z <- if (methods::is(object@randomizermod, "randsdvine")) {
      randsdvine_sample(V1, V2, object@basecopula, object@randomizermod)
    } else {
      randmixture_sample(V1, V2, object@randomizermod)
    }
    U1 <- udpsi(object@udp1, V1, Z[, "Z1"])
    U2 <- udpsi(object@udp2, V2, Z[, "Z2"])
  }
  cbind(U1 = U1, U2 = U2)
}

# Column of udpinverse(x, v, prob = TRUE)'s (sorted-ascending) pre-image
# matrix M that matches u, within a numerical tolerance. u is assumed
# consistent with v (v == udptrans(x, u) for the same x), so a match should
# always exist -- the tolerance only absorbs udpinverse()'s own imprecision,
# not genuine ambiguity. It matches the ~1e-4 accuracy of the spline-based
# classes (udplegendre, udplegendrebex, udpcosinebex), whose pre-images go
# through an interpolated F^{-1}: the round trip u -> v -> pre-images is
# usually good to 1e-9, but near a turning point of g, where roots move
# like the square root of an error in v, it can be off by 1e-6 or more.
match_preimage <- function(M, u, tol = 1e-4) {
  d <- abs(M - u)
  d[is.na(d)] <- Inf
  j <- max.col(-d, ties.method = "first")
  mismatch <- d[cbind(seq_along(u), j)]
  if (any(mismatch > tol)) {
    stop(
      "'u1'/'u2' is not consistent with 'udp1'/'udp2' (no matching pre-image found).",
      call. = FALSE
    )
  }
  j
}

# The selection-probability interval [a, b] in [0, 1] (Z-space) that
# udpsi() assigns to column j of a udpinverse(..., prob = TRUE) result: the
# running sum of that row's selection probabilities up to (a, exclusive) and
# including (b) column j. Mirrors the cumulative construction udpsi() uses
# internally to turn a randomizer draw into a column choice. The ends are
# clamped into [0, 1]: a running sum of probabilities that add to 1 can
# overshoot to 1 + 2e-16 by rounding, and rvinecopulib rejects arguments
# outside [0, 1].
preimage_interval <- function(P, j) {
  P[is.na(P)] <- 0
  k <- ncol(P)
  n <- nrow(P)
  cum <- pmin(P %*% upper.tri(matrix(0, k, k), diag = TRUE), 1)
  b <- cum[cbind(seq_len(n), j)]
  a <- ifelse(j == 1L, 0, cum[cbind(seq_len(n), pmax(j - 1L, 1L))])
  list(a = a, b = b)
}

# CDF of a copula (parCopula or bicop_dist), exact on the boundary of the
# unit square: C(a, 0) = C(0, b) = 0, C(a, 1) = a, C(1, b) = b. rvinecopulib
# clips its arguments to [1e-10, 1 - 1e-10], so its own values there are off
# by ~1e-10 -- enough to stop rectangle probabilities over a partition of
# [0, 1]^2 summing to exactly 1.
exact_copula_cdf <- function(a, b, cop) {
  out <- numeric(length(a))
  # known on the boundary, so call the copula only at interior points
  hi <- a >= 1 & b > 0
  out[hi] <- b[hi]
  hi <- b >= 1 & a > 0
  out[hi] <- a[hi]
  i <- which(a > 0 & a < 1 & b > 0 & b < 1)
  if (length(i)) {
    out[i] <- basecopula_cdf(cop, a[i], b[i])
  }
  out
}

# h-function P(Z <= z | E = e) of a bicop_dist with variables (E, Z) if
# z_second is TRUE, (Z, E) otherwise, exact at z = 0 and z = 1 for the same
# reason as exact_copula_cdf().
exact_hfunc <- function(e, z, cop, z_second = TRUE) {
  h <- if (z_second) {
    rvinecopulib::hbicop(cbind(e, z), cond_var = 1, family = cop)
  } else {
    rvinecopulib::hbicop(cbind(z, e), cond_var = 2, family = cop)
  }
  h[z <= 0] <- 0
  h[z >= 1] <- 1
  h
}

# The joint conditional CDF of the randomizers, F(z1, z2 | V1 = v1, V2 = v2),
# as a function F(z1, z2, idx) evaluating it at the (v1, v2) pairs indexed by
# idx (all of them, in order, by default). Everything that depends only on
# (v1, v2) is computed once here, however many times F is then called.
#   NULL        independent uniform randomizers: z1 * z2.
#   randsdvine  the tree-3 copula applied to the two tree-2 h-functions, each
#               conditioned on the matching tree-1 h-function of (v1, v2):
#               e21 = P(V2 <= v2 | V1 = v1), e12 = P(V1 <= v1 | V2 = v2), the
#               same quantities randsdvine_sample() uses for simulation.
#               copZ1V2_V1 has variables (Z1, V2), so its h-function for Z1
#               conditions on the second one; copV1Z2_V2 has (V1, Z2), so its
#               h-function for Z2 conditions on the first.
#   randmixture selector(v1, v2) picks cop1 or cop2, and (Z1, Z2) given that
#               choice is an ordinary draw from the picked copula, so F is
#               that copula's own CDF.
randomizer_cdf <- function(v1, v2, basecopula, randomizermod) {
  all_idx <- seq_along(v1)
  if (is.null(randomizermod)) {
    return(function(z1, z2, idx = all_idx) z1 * z2)
  }
  if (methods::is(randomizermod, "randsdvine")) {
    e21 <- basecopula_cond_cdf(basecopula, v1, v2, 1L)
    e12 <- basecopula_cond_cdf(basecopula, v2, v1, 2L)
    return(function(z1, z2, idx = all_idx) {
      exact_copula_cdf(
        exact_hfunc(e21[idx], z1, randomizermod@copZ1V2_V1, z_second = FALSE),
        exact_hfunc(e12[idx], z2, randomizermod@copV1Z2_V2, z_second = TRUE),
        randomizermod@copZ1Z2_V1V2
      )
    })
  }
  sel <- randomizermod@selector(v1, v2)
  validate_selector_result(sel, length(v1))
  function(z1, z2, idx = all_idx) {
    s <- sel[idx]
    out <- numeric(length(z1))
    if (any(s)) {
      out[s] <- exact_copula_cdf(z1[s], z2[s], randomizermod@cop1)
    }
    if (any(!s)) {
      out[!s] <- exact_copula_cdf(z1[!s], z2[!s], randomizermod@cop2)
    }
    out
  }
}

# w(u1, u2): the ratio of (1) the probability that the joint stochastic
# inverse (udpsi(udp1, V1, Z1), udpsi(udp2, V2, Z2)) lands in the cell
# containing (u1, u2), given V1 = v1, V2 = v2, under randomizermod, to (2)
# the same probability under independent randomizers -- a rectangle
# probability under the joint conditional law of (Z1, Z2) | (V1, V2)
# (numerator, via randomizer_cdf()) over the product of the two
# marginal selection probabilities already in udpinverse()'s "prob"
# attribute (denominator). See McNeil and Nešlehová (2026), arXiv:2607.07174,
# also cited via dbsicopula()'s @references.
bsicopula_weight <- function(u1, u2, v1, v2, udp1, udp2, basecopula, randomizermod) {
  M1 <- udpinverse(udp1, v1, prob = TRUE)
  M2 <- udpinverse(udp2, v2, prob = TRUE)
  P1 <- attr(M1, "prob")
  P2 <- attr(M2, "prob")

  j1 <- match_preimage(M1, u1)
  j2 <- match_preimage(M2, u2)
  I1 <- preimage_interval(P1, j1)
  I2 <- preimage_interval(P2, j2)

  denom <- P1[cbind(seq_along(u1), j1)] * P2[cbind(seq_along(u2), j2)]

  joint_cdf <- randomizer_cdf(v1, v2, basecopula, randomizermod)

  numer <- joint_cdf(I1$b, I2$b) - joint_cdf(I1$a, I2$b) -
    joint_cdf(I1$b, I2$a) + joint_cdf(I1$a, I2$a)

  numer / denom
}

#' Density of a bivariate stochastic inversion copula
#'
#' The population density of a \linkS4class{bsicopula}. When `randomizermod`
#' is `NULL` (independent stochastic inversion), this is simply
#' `c_V(udptrans(udp1, u1), udptrans(udp2, u2))`, `c_V` being the density of
#' `basecopula`: `udpsi()` always returns an exact pre-image, and the
#' `1 / |T'|` Jacobian of the change of variables cancels exactly against
#' the `1 / |T'|` pre-image selection probability, since a
#' \linkS4class{udp} transformation's pre-image selection probabilities sum
#' to `1` at every `v` by construction.
#'
#' When `randomizermod` is a \linkS4class{randsdvine} or a
#' \linkS4class{randmixture}, that cancellation is no longer exact -- `Z1` now
#' depends on `V2` (and vice versa), not just on its own `V_i` -- and an
#' extra multiplicative weight `w(u1, u2)` corrects for it: the ratio of
#' the probability that the joint stochastic inverse lands in the same
#' pair of pre-image branches as `(u1, u2)`, given `(V1, V2)`, under
#' `randomizermod`, to the same probability under independent randomizers.
#' `w` is `1` identically when `randomizermod` is `NULL`.
#'
#' `udp1` and `udp2` are each piecewise smooth with finitely many
#' breakpoints; at a breakpoint of either, the number of
#' pre-image branches drops (two branches meet at the same point), so the
#' density has a genuine jump discontinuity there -- harmless for a density
#' (a finite set of points has measure zero) but worth knowing about if
#' evaluating `dbsicopula()` on a grid that happens to include one:
#' `udptrans()`/`udpinverse()` are recomputed fresh from `u1`, `u2` at every
#' call, so the value returned at a breakpoint is the (well-defined) value
#' for the merged branch there, not an arbitrary one-sided limit.
#'
#' @param u1,u2 numeric vectors with values in `[0, 1]`, of equal length or
#'   with a length-1 argument recycled to the length of the other. `u1` may
#'   instead be a two-column matrix (or data frame) with `u2` omitted, as for
#'   [rvinecopulib::dbicop()]: `dbsicopula(U, object = bc)`, or
#'   `dbsicopula(U, bc)`.
#' @param object an object of class \linkS4class{bsicopula}.
#'
#' @return A numeric vector of density values, one per observation.
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), udpcosine(2), udpcosine(3))
#'   dbsicopula(c(0.2, 0.5), c(0.3, 0.5), bc)
#' }
dbsicopula <- function(u1, u2 = NULL, object) {
  a <- bsicopula_args(u1, u2, object, if (missing(object)) NULL else object)
  dbsicopula_eval(a$u1, a$u2, a$object)
}

# Common argument handling of dbsicopula(), pbsicopula() and hbsicopula():
# (u1, u2) as two vectors, or u1 as a two-column matrix with u2 omitted (and
# then, when 'object' is omitted too, the object may sit in the u2 position:
# dbsicopula(U, bc)); a length-1 vector is recycled to the length of the
# other. 'range' additionally requires every value to be a number in [0, 1].
bsicopula_args <- function(u1, u2, object, object_given, range = FALSE) {
  if (is.null(object_given)) {
    if (!methods::is(u2, "bsicopula")) {
      stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
    }
    object <- u2
    u2 <- NULL
  }
  if (!methods::is(object, "bsicopula")) {
    stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
  }
  if (is.null(u2)) {
    if (!(is.matrix(u1) || is.data.frame(u1)) || ncol(u1) != 2L) {
      stop("when 'u2' is omitted, 'u1' must be a two-column matrix.", call. = FALSE)
    }
    u1 <- as.matrix(u1)
    u2 <- as.numeric(u1[, 2L])
    u1 <- as.numeric(u1[, 1L])
  } else {
    u1 <- as.numeric(u1)
    u2 <- as.numeric(u2)
    n <- max(length(u1), length(u2))
    if (length(u1) == 1L && n != 1L) u1 <- rep(u1, n)
    if (length(u2) == 1L && n != 1L) u2 <- rep(u2, n)
    if (length(u1) != length(u2)) {
      stop("'u1' and 'u2' must have the same length (or one of them length 1).", call. = FALSE)
    }
  }
  if (range && (anyNA(u1) || anyNA(u2) || any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1))) {
    stop("'u1' and 'u2' must be numbers in [0, 1], without missing values.", call. = FALSE)
  }
  list(u1 = u1, u2 = u2, object = object)
}

# Clamp carrier values V into [f, 1 - f], but only for observations whose u
# is itself inside [f, 1 - f]. An observation that close to 0 or 1 sits
# next to one of the fixed endpoints every udp maps to 0 or 1 whatever its
# parameters, so its extreme V is genuine data, not the product of a moving
# breakpoint, and is left exact. For pseudo-observations (all u in
# [1/(n+1), n/(n+1)]) and f = 1/(2n), every observation is clamped.
clamp_v <- function(V, u, f) {
  inside <- u >= f & u <= 1 - f
  V[inside] <- pmin(pmax(V[inside], f), 1 - f)
  V
}

# dbsicopula()'s computation, minus argument checking, with one extra option
# for fitbsicopula(): when vfloor is non-NULL, the carrier values are clamped
# into [vfloor, 1 - vfloor] (see clamp_v()) before the base copula density is
# evaluated -- and only there. The weight w(u1, u2) still sees the exact
# V1, V2, since match_preimage() needs u_i and V_i to be consistent.
dbsicopula_eval <- function(u1, u2, object, vfloor = NULL) {
  # rounding can put a v-transform a few ulps outside [0, 1]
  V1 <- pmin(pmax(udptrans(object@udp1, u1), 0), 1)
  V2 <- pmin(pmax(udptrans(object@udp2, u2), 0), 1)
  cV <- if (is.null(vfloor)) {
    basecopula_density(object@basecopula, V1, V2)
  } else {
    basecopula_density(
      object@basecopula, clamp_v(V1, u1, vfloor), clamp_v(V2, u2, vfloor)
    )
  }

  if (is.null(object@randomizermod)) {
    return(cV)
  }
  w <- bsicopula_weight(
    u1, u2, V1, V2, object@udp1, object@udp2, object@basecopula, object@randomizermod
  )
  cV * w
}

#' Plot method for the bsicopula class
#'
#' Draws the bivariate density of a \linkS4class{bsicopula} object
#' ([dbsicopula()]) over a regular grid, as a contour plot (`type =
#' "contour"`, the default) or a perspective plot (`type = "persp"`).
#'
#' Densities built from tail-dependent copulas or many-branch `udp1`/`udp2`
#' transforms are often extremely right-skewed (a small high-density region
#' against a large near-zero one). `type = "contour"` handles this in two
#' ways: contour `levels` default to `quantile(dens, probs)` -- by default
#' the median and up -- rather than [graphics::contour()]'s own evenly
#' spaced `pretty()` levels, which would fall almost entirely above the
#' bulk of the mass and leave the plot looking empty; and the filled
#' background is coloured by each grid point's own rank among all evaluated
#' density values (`rank(dens) / length(dens)`), not its raw value, so the
#' colour scale is spread evenly across the whole plot instead of collapsing
#' into a few small high-density patches against an undifferentiated
#' background. The trade-off is that the fill colour is ordinal, not a
#' literal density scale -- the contour lines carry the actual values.
#' Before ranking, every value below the lowest contour `level` is floored
#' to that level, so they all tie for the same (lowest) rank: without this,
#' near-degenerate copulas (correlation close to `-1` or `1`, say) can leave
#' a speckle of tiny but numerically nonzero density values scattered across
#' what should be a uniform background, which `rank()` -- having no notion
#' of "negligible" -- would otherwise render as visible texture.
#'
#' @param x an object of class \linkS4class{bsicopula}.
#' @param type `"contour"` (the default) or `"persp"`.
#' @param n number of grid points per axis. The grid uses cell midpoints
#'   `(i - 0.5) / n`, which generically avoids landing exactly on a
#'   breakpoint of `udp1` or `udp2` -- see [dbsicopula()] for why the
#'   density has a jump discontinuity there.
#' @param probs for `type = "contour"`: probabilities passed to
#'   [stats::quantile()] on the evaluated density to give the default
#'   contour `levels`. Ignored if `levels` is supplied directly.
#' @param levels for `type = "contour"`: contour levels, overriding the
#'   `probs`-based default.
#' @param col for `type = "contour"`: the fill colour palette, recycled
#'   across the rank-transformed density (see Details); passed to
#'   [graphics::image()] as `col`.
#' @param drawlabels for `type = "contour"`: whether to label the contour
#'   lines with their level. Default `FALSE`, since the `probs`-based
#'   levels are typically numerous enough that labels overlap and clutter
#'   the plot.
#' @param xlab,ylab axis labels.
#' @param ... further graphical parameters passed to [graphics::image()]
#'   (`type = "contour"`) or [graphics::persp()] (`type = "persp"`).
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), udpcosine(2), udpcosine(3))
#'   plot(bc)
#'   plot(bc, type = "persp")
#' }
setMethod("plot", c(x = "bsicopula", y = "missing"),
  function(x, type = c("contour", "persp"), n = 100, xlab = "u1", ylab = "u2",
           probs = c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999), levels = NULL,
           col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE), drawlabels = FALSE, ...) {
    type <- match.arg(type)
    grid <- (seq_len(n) - 0.5) / n
    U <- expand.grid(u1 = grid, u2 = grid)
    dens <- matrix(dbsicopula(U$u1, U$u2, x), n, n)
    if (type == "contour") {
      if (is.null(levels)) {
        levels <- sort(unique(stats::quantile(dens, probs)))
      }
      # Rank-transform for the fill colour (see Details), but first floor
      # every value below the lowest drawn contour level to that level, so
      # they all tie for the same (lowest) rank/colour. Without this,
      # everything below the lowest level is nominally uniform "background"
      # but the raw density there often isn't identically 0 -- near-
      # degenerate copulas in particular (e.g. correlation close to +-1)
      # leave a speckle of tiny, numerically noisy nonzero values -- and
      # rank() has no notion of "negligible", so those specks get
      # perceptibly different colours from an otherwise flat background.
      dens_for_rank <- pmax(dens, min(levels))
      rankdens <- matrix(rank(dens_for_rank) / length(dens_for_rank), n, n)
      dots <- list(...)
      if (is.null(dots$asp)) dots$asp <- 1
      if (is.null(dots$xaxs)) dots$xaxs <- "i"
      if (is.null(dots$yaxs)) dots$yaxs <- "i"
      # pty = "s" makes the plotting region square in physical inches;
      # without it, asp = 1 stretches one axis's displayed range past its
      # data limits on any device or panel that isn't already exactly square.
      op <- graphics::par(pty = "s")
      on.exit(graphics::par(op))
      do.call(image, c(
        list(x = grid, y = grid, z = rankdens, col = col, xlab = xlab, ylab = ylab), dots
      ))
      contour(grid, grid, dens, levels = levels, drawlabels = drawlabels, add = TRUE, col = "grey20")
    } else {
      persp(grid, grid, dens, xlab = xlab, ylab = ylab, zlab = "density", ...)
    }
  }
)
