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
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
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
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
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

# Density of a base copula that is either a parCopula or a bicop_dist
# object, dispatching to the matching package's own density function.
basecopula_density <- function(v1, v2, basecopula) {
  if (is_parCopula(basecopula)) {
    copula::dCopula(cbind(v1, v2), basecopula)
  } else {
    rvinecopulib::dbicop(cbind(v1, v2), basecopula)
  }
}

# Column of udpinverse(x, v, prob = TRUE)'s (sorted-ascending) pre-image
# matrix M that matches u, within a numerical tolerance. u is assumed
# consistent with v (v == udptrans(x, u) for the same x), so a match should
# always exist -- the tolerance only absorbs udpinverse()'s own root-finding
# imprecision, not genuine ambiguity.
match_preimage <- function(M, u, tol = 1e-6) {
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
# internally to turn a randomizer draw into a column choice.
preimage_interval <- function(P, j) {
  P[is.na(P)] <- 0
  k <- ncol(P)
  n <- nrow(P)
  cum <- P %*% upper.tri(matrix(0, k, k), diag = TRUE)
  b <- cum[cbind(seq_len(n), j)]
  a <- ifelse(j == 1L, 0, cum[cbind(seq_len(n), pmax(j - 1L, 1L))])
  list(a = a, b = b)
}

# Joint conditional CDF F(z1, z2 | V1 = v1, V2 = v2) implied by the
# simplified D-vine: the tree-3 copula applied to the two tree-2
# h-functions, each conditioned on the matching tree-1 h-function of
# (v1, v2). e21 = P(V2 <= v2 | V1 = v1), e12 = P(V1 <= v1 | V2 = v2) -- the
# same quantities sdvine_sample() computes for simulation.
sdvine_cdf <- function(z1, z2, e21, e12, randomizermod) {
  a <- rvinecopulib::hbicop(cbind(e21, z1), cond_var = 1, family = randomizermod@copZ1V2_V1)
  b <- rvinecopulib::hbicop(cbind(e12, z2), cond_var = 1, family = randomizermod@copV1Z2_V2)
  rvinecopulib::pbicop(cbind(a, b), family = randomizermod@copZ1Z2_V1V2)
}

# w(u1, u2): the ratio of (1) the probability that the joint stochastic
# inverse (udpsi(udp1, V1, Z1), udpsi(udp2, V2, Z2)) lands in the cell
# containing (u1, u2), given V1 = v1, V2 = v2, under randomizermod, to (2)
# the same probability under independent randomizers -- a rectangle
# probability under the joint conditional law of (Z1, Z2) | (V1, V2)
# (numerator, via sdvine_cdf()'s closed form) over the product of the two
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

  e21 <- rvinecopulib::hbicop(cbind(v1, v2), cond_var = 1, family = basecopula)
  e12 <- rvinecopulib::hbicop(cbind(v2, v1), cond_var = 1, family = basecopula)

  numer <- sdvine_cdf(I1$b, I2$b, e21, e12, randomizermod) -
    sdvine_cdf(I1$a, I2$b, e21, e12, randomizermod) -
    sdvine_cdf(I1$b, I2$a, e21, e12, randomizermod) +
    sdvine_cdf(I1$a, I2$a, e21, e12, randomizermod)

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
#' When `randomizermod` is an \linkS4class{sdvine}, that cancellation is no
#' longer exact -- `Z1` now depends on `V2` (and vice versa) through the
#' D-vine, not just on its own `V_i` -- and an extra multiplicative weight
#' `w(u1, u2)` corrects for it: the ratio of the probability that the joint
#' stochastic inverse lands in the same pair of pre-image branches as
#' `(u1, u2)`, given `(V1, V2)`, under `randomizermod`, to the same
#' probability under independent randomizers. `w` is `1` identically when
#' `randomizermod` is `NULL`.
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
#' @param u1,u2 numeric vectors of equal length with values in `[0, 1]`.
#' @param object an object of class \linkS4class{bsicopula}.
#'
#' @return A numeric vector, the same length as `u1`, of density values.
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
dbsicopula <- function(u1, u2, object) {
  if (!methods::is(object, "bsicopula")) {
    stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
  }
  u1 <- as.numeric(u1)
  u2 <- as.numeric(u2)
  if (length(u1) != length(u2)) {
    stop("'u1' and 'u2' must have the same length.", call. = FALSE)
  }

  V1 <- udptrans(object@udp1, u1)
  V2 <- udptrans(object@udp2, u2)
  cV <- basecopula_density(V1, V2, object@basecopula)

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
#' @param x an object of class \linkS4class{bsicopula}.
#' @param type `"contour"` (the default) or `"persp"`.
#' @param n number of grid points per axis. The grid uses cell midpoints
#'   `(i - 0.5) / n`, which generically avoids landing exactly on a
#'   breakpoint of `udp1` or `udp2` -- see [dbsicopula()] for why the
#'   density has a jump discontinuity there.
#' @param xlab,ylab axis labels.
#' @param ... further graphical parameters passed to [graphics::contour()]
#'   (`type = "contour"`) or [graphics::persp()] (`type = "persp"`), such as
#'   `levels` to override the default contour levels. Densities built from
#'   tail-dependent copulas or many-branch `udp1`/`udp2` transforms are
#'   often extremely right-skewed (a small high-density region against a
#'   large near-zero one), so `type = "contour"` defaults to the deciles of
#'   the evaluated density (`quantile(dens, seq(0.1, 0.9, by = 0.1))`)
#'   rather than [graphics::contour()]'s own evenly spaced `pretty()`
#'   levels, which tend to fall almost entirely above the bulk of the mass
#'   in that case and leave the plot looking empty.
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
  function(x, type = c("contour", "persp"), n = 100, xlab = "u1", ylab = "u2", ...) {
    type <- match.arg(type)
    grid <- (seq_len(n) - 0.5) / n
    U <- expand.grid(u1 = grid, u2 = grid)
    dens <- matrix(dbsicopula(U$u1, U$u2, x), n, n)
    if (type == "contour") {
      dots <- list(...)
      if (is.null(dots$levels)) {
        dots$levels <- sort(unique(stats::quantile(dens, probs = seq(0.1, 0.9, by = 0.1))))
      }
      do.call(contour, c(list(x = grid, y = grid, z = dens, xlab = xlab, ylab = ylab), dots))
    } else {
      persp(grid, grid, dens, xlab = xlab, ylab = ylab, zlab = "density", ...)
    }
  }
)
