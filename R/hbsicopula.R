# CDF and h-functions of a bivariate stochastic inversion copula.
#
# With independent randomizers (randomizermod = NULL) the copula of
# (U1, U2) has density c_V(T1(u1), T2(u2)), c_V the density of the base
# copula. If T1 and T2 are piecewise linear, every integral is then a finite
# sum: piece k of T maps u in [a_k, b_k] linearly onto a v-interval [c_k, d_k]
# and is selected, given V = v in [c_k, d_k], with probability w_k (see
# udplinpieces()), so
#
#   C(u1, u2)     = sum_{j,k} w1_j w2_k * mass of C_V on I1_j(u1) x I2_k(u2)
#   h(u2 | u1)    = sum_k w2_k * [H(hi_k | v1) - H(lo_k | v1)],  v1 = T1(u1)
#
# where I_k(u) = [lo_k, hi_k] is the v-interval the piece sweeps as its
# argument runs up to u and H(x | v1) = P(V2 <= x | V1 = v1) is the base
# copula's h-function. Only the margin integrated over need be piecewise
# linear for h; the conditioning margin enters through v1 = T1(u1) alone.
# For a linear v-transform on both margins this is Proposition S3 of the
# supplement to Dias, Han and McNeil, "GARCH copulas, v-transforms and
# D-vines for stochastic volatility".

# P(X <= x | G = given) under the base copula, X being the carrier other than
# carrier 'given_var' (1 or 2), exact at x = 0 and x = 1 and clipped to
# [0, 1]; the interior values come from basecopula_h().
basecopula_cond_cdf <- function(cop, given, x, given_var) {
  out <- numeric(length(x))
  out[x >= 1] <- 1
  i <- which(x > 0 & x < 1) # at x = 0 and 1 the value is known, so skip the call
  if (length(i)) {
    out[i] <- basecopula_h(cop, given[i], x[i], given_var)
  }
  pmin(pmax(out, 0), 1)
}

# The inverse in x of basecopula_cond_cdf(): the x with P(X <= x | G = given)
# equal to q, from basecopula_hinv().
basecopula_cond_quantile <- function(cop, given, q, given_var) {
  basecopula_hinv(cop, given, pmin(pmax(q, 0), 1), given_var)
}

# The v-interval swept by piece k of P as the argument runs from 0 up to u,
# v = T(u): empty before the piece, the whole [c_k, d_k] after it, and
# between, the image of [a_k, u].
swept_interval <- function(P, k, u, v) {
  incr <- P$incr[k]
  lo <- if (incr) rep(P$c[k], length(u)) else pmin(pmax(v, P$c[k]), P$d[k])
  hi <- if (incr) pmin(pmax(v, P$c[k]), P$d[k]) else rep(P$d[k], length(u))
  above <- u >= P$b[k]
  below <- u <= P$a[k]
  lo[above] <- P$c[k]
  hi[above] <- P$d[k]
  lo[below] <- P$c[k]
  hi[below] <- P$c[k]
  list(lo = lo, hi = hi)
}

# Exact CDF for piecewise-linear margins: the sum over pairs of pieces of the
# product of their weights and the base copula's mass of the rectangle of
# swept v-intervals.
pbsicopula_exact <- function(u1, u2, object) {
  P1 <- udplinpieces(object@udp1)
  P2 <- udplinpieces(object@udp2)
  v1 <- udptrans(object@udp1, u1)
  v2 <- udptrans(object@udp2, u2)
  cop <- object@basecopula
  out <- numeric(length(u1))
  for (j in seq_len(P1$K)) {
    A <- swept_interval(P1, j, u1, v1)
    for (k in seq_len(P2$K)) {
      B <- swept_interval(P2, k, u2, v2)
      i <- which(A$hi > A$lo & B$hi > B$lo)
      if (!length(i)) next
      rect <- exact_copula_cdf(A$hi[i], B$hi[i], cop) - exact_copula_cdf(A$lo[i], B$hi[i], cop) -
        exact_copula_cdf(A$hi[i], B$lo[i], cop) + exact_copula_cdf(A$lo[i], B$lo[i], cop)
      out[i] <- out[i] + P1$w[j] * P2$w[k] * rect
    }
  }
  pmin(pmax(out, 0), 1)
}

# Exact h-function P(U_t <= u_t | U_g = u_g) for a piecewise-linear target
# margin t (pieces Pt), the conditioning margin g entering through its
# carrier value v_g = T_g(u_g).
h_exact <- function(u_t, v_t, v_g, Pt, cop, given_var) {
  out <- numeric(length(u_t))
  for (k in seq_len(Pt$K)) {
    B <- swept_interval(Pt, k, u_t, v_t)
    i <- which(B$hi > B$lo)
    if (!length(i)) next
    out[i] <- out[i] + Pt$w[k] * (basecopula_cond_cdf(cop, v_g[i], B$hi[i], given_var) -
      basecopula_cond_cdf(cop, v_g[i], B$lo[i], given_var))
  }
  pmin(pmax(out, 0), 1)
}

# Exact inverse of h_exact() in u_t. The mass of piece k given v_g is
# w_k * [H(d_k | v_g) - H(c_k | v_g)] (just w_k when the piece maps onto
# [0, 1]); the level p falls in the first piece whose cumulative mass reaches
# it, and within that piece solves the linear equation for the swept
# interval, so one inverse base-copula h-function finishes the job.
h_inverse_exact <- function(p, v_g, Pt, cop, given_var) {
  n <- length(p)
  if (!n) {
    return(numeric(0))
  }
  K <- Pt$K
  Hat <- function(x) {
    if (x <= 0) rep(0, n) else if (x >= 1) rep(1, n) else basecopula_cond_cdf(cop, v_g, rep(x, n), given_var)
  }
  Hc <- Hd <- cm <- matrix(0, n, K)
  for (k in seq_len(K)) {
    Hc[, k] <- Hat(Pt$c[k])
    Hd[, k] <- Hat(Pt$d[k])
    cm[, k] <- Pt$w[k] * (Hd[, k] - Hc[, k]) + if (k > 1L) cm[, k - 1L] else 0
  }
  piece <- if (K > 1L) 1L + rowSums(p > cm[, -K, drop = FALSE]) else rep(1L, n)
  at <- cbind(seq_len(n), piece)
  before <- ifelse(piece == 1L, 0, cm[cbind(seq_len(n), pmax(piece - 1L, 1L))])
  r <- p - before
  w <- Pt$w[piece]
  incr <- Pt$incr[piece]
  level <- ifelse(incr, Hc[at] + r / w, Hd[at] - r / w)
  v <- basecopula_cond_quantile(cop, v_g, level, given_var)
  u <- ifelse(incr, Pt$a[piece] + w * (v - Pt$c[piece]), Pt$a[piece] + w * (Pt$d[piece] - v))
  pmin(pmax(u, Pt$a[piece]), Pt$b[piece])
}

# ---- Quadrature tier: margins that are not piecewise linear ----------------
#
# With independent randomizers the density is c_V(T1(u1), T2(u2)), so
#
#   h(u_t | v_g) = integral_0^{u_t} c_V(v_g, T_t(x)) dx
#
# is a sum of integrals over the cells of udpbreaks(T_t), on each of which the
# integrand is smooth up to its ends (where it can be singular: T' may be
# infinite there). They are done by tanh-sinh quadrature, which is insensitive
# to endpoint singularities. The CDF is the integral over one margin of the
# h-function of the other, conditioning on that margin's carrier value; the
# inner margin may be piecewise linear (exact h) or not.

# Tanh-sinh nodes on (0, 1): x with the complement 1 - x computed separately,
# so nodes near 1 keep their relative accuracy, and weights summing to ~1.
tanh_sinh_nodes <- function(N) {
  t <- seq(-3.2, 3.2, length.out = N)
  s <- pi * sinh(t)
  list(
    x = 1 / (1 + exp(-s)), xc = 1 / (1 + exp(s)),
    w = (t[2] - t[1]) * pi * cosh(t) / ((1 + exp(-s)) * (1 + exp(s)))
  )
}

# The integral of f over [a, upper[i]] for each i. f(t, idx) takes the nodes t
# and the index idx of the row each belongs to, and returns the integrand
# there. Rows go through in chunks to bound memory.
cell_quad <- function(a, upper, N, f, chunk = 2e5) {
  n <- length(upper)
  a <- rep_len(a, n)
  out <- numeric(n)
  live <- which(upper > a)
  if (!length(live)) {
    return(out)
  }
  nd <- tanh_sinh_nodes(N)
  per <- max(1L, floor(chunk / N))
  for (start in seq(1L, length(live), by = per)) {
    i <- live[start:min(start + per - 1L, length(live))]
    m <- length(i)
    L <- matrix(upper[i] - a[i], N, m, byrow = TRUE)
    X <- matrix(nd$x, N, m)
    XC <- matrix(nd$xc, N, m)
    t <- ifelse(X < 0.5, matrix(a[i], N, m, byrow = TRUE) + L * X, matrix(upper[i], N, m, byrow = TRUE) - L * XC)
    val <- f(as.vector(t), rep(i, each = N))
    val[!is.finite(val)] <- 0
    out[i] <- L[1L, ] * colSums(matrix(val, N) * nd$w)
  }
  out
}

# The integrand of the h-function: the density of the target margin at t given
# the conditioning margin's value g[idx]. With independent randomizers this is
# c_V(v_g, T_t(t)), g holding the carrier values v_g. With a randomizer model
# the density carries the weight w(u1, u2), which depends on the branches
# selected, so g holds the conditioning value u_g and the integrand is
# dbsicopula() itself.
h_integrand <- function(g, udp_t, object, given_var) {
  if (is.null(object@randomizermod)) {
    cop <- object@basecopula
    return(function(t, idx) {
      vt <- pmin(pmax(udptrans(udp_t, t), 0), 1)
      if (given_var == 1L) basecopula_density(cop, g[idx], vt) else basecopula_density(cop, vt, g[idx])
    })
  }
  function(t, idx) {
    if (given_var == 1L) dbsicopula_eval(g[idx], t, object) else dbsicopula_eval(t, g[idx], object)
  }
}

# h(u_t | g) = integral_0^{u_t} of the density over the cells of udp_t; see
# h_integrand() for what g is.
h_quad <- function(u_t, g, udp_t, object, given_var, N) {
  br <- udpbreaks(udp_t)
  dens <- h_integrand(g, udp_t, object, given_var)
  out <- numeric(length(u_t))
  for (k in seq_len(length(br) - 1L)) {
    out <- out + cell_quad(br[k], pmin(pmax(u_t, br[k]), br[k + 1L]), N, dens)
  }
  pmin(pmax(out, 0), 1)
}

# Inverse of h_quad() in u_t. Cumulative cell masses locate the cell that
# holds the level; within it a bracketed Newton iteration (the derivative is
# the density, bisection when a step leaves the bracket) solves for u_t.
h_inverse_quad <- function(p, g, udp_t, object, given_var, N, tol = 1e-12, maxit = 100L) {
  n <- length(p)
  if (!n) {
    return(numeric(0))
  }
  br <- udpbreaks(udp_t)
  K <- length(br) - 1L
  p <- pmin(pmax(p, 0), 1)
  dens <- h_integrand(g, udp_t, object, given_var)
  cm <- matrix(0, n, K)
  for (k in seq_len(K)) {
    cm[, k] <- cell_quad(br[k], rep(br[k + 1L], n), N, dens) + if (k > 1L) cm[, k - 1L] else 0
  }
  piece <- if (K > 1L) pmin(1L + rowSums(p > cm[, -K, drop = FALSE]), K) else rep(1L, n)
  before <- ifelse(piece == 1L, 0, cm[cbind(seq_len(n), pmax(piece - 1L, 1L))])
  lo <- br[piece]
  hi <- br[piece + 1L]
  mass <- cm[cbind(seq_len(n), piece)] - before
  u <- lo + (hi - lo) * pmin(pmax((p - before) / pmax(mass, 1e-300), 0), 1)
  active <- which(p > 0 & p < 1)
  for (it in seq_len(maxit)) {
    if (!length(active)) break
    r <- active
    f <- before[r] + cell_quad(br[piece[r]], u[r], N, function(t, idx) dens(t, r[idx])) - p[r]
    below <- f < 0
    lo[r[below]] <- u[r[below]]
    hi[r[!below]] <- u[r[!below]]
    done <- abs(f) < tol | (hi[r] - lo[r]) < 1e-15
    step <- u[r] - f / pmax(dens(u[r], r), 1e-300)
    unew <- ifelse(step > lo[r] & step < hi[r], step, (lo[r] + hi[r]) / 2)
    u[r] <- ifelse(done, u[r], unew)
    active <- r[!done]
  }
  pmin(pmax(u, 0), 1)
}

# h(u_t | g) for the target margin t of 'object', exact if the target is
# piecewise linear and the randomizers independent, by quadrature otherwise.
# g is the conditioning margin's carrier value v_g in the first case, its
# argument u_g in the second (see h_integrand()).
h_any <- function(u_t, v_t, g, udp_t, object, given_var, N) {
  Pt <- udplinpieces(udp_t)
  if (is.null(Pt) || !is.null(object@randomizermod)) {
    h_quad(u_t, g, udp_t, object, given_var, N)
  } else {
    h_exact(u_t, v_t, g, Pt, object@basecopula, given_var)
  }
}

# CDF when a margin is not piecewise linear, or the randomizers are not
# independent: integrate over one (outer) margin the h-function of the other,
# conditioning on the outer margin's carrier value, or, with a randomizer,
# on its argument. Two quadratures nest when both margins need one: N^2 nodes
# per pair of cells.
pbsicopula_quad <- function(u1, u2, object, N) {
  o <- if (is.null(udplinpieces(object@udp2)) || !is.null(object@randomizermod)) 2L else 1L
  udp_o <- slot(object, paste0("udp", o))
  udp_i <- slot(object, paste0("udp", 3L - o))
  u_o <- if (o == 1L) u1 else u2
  u_i <- if (o == 1L) u2 else u1
  v_i <- udptrans(udp_i, u_i)
  rand <- !is.null(object@randomizermod)
  inner <- function(t, idx) {
    h_any(u_i[idx], v_i[idx], if (rand) t else udptrans(udp_o, t), udp_i, object, o, N)
  }
  br <- udpbreaks(udp_o)
  out <- numeric(length(u_o))
  for (k in seq_len(length(br) - 1L)) {
    out <- out + cell_quad(br[k], pmin(pmax(u_o, br[k]), br[k + 1L]), N, inner)
  }
  pmin(pmax(out, 0), 1)
}

check_nodes <- function(nodes) {
  if (!is.numeric(nodes) || length(nodes) != 1L || is.na(nodes) || nodes < 3 || nodes != round(nodes)) {
    stop("'nodes' must be a single integer of at least 3.", call. = FALSE)
  }
}

#' CDF of a bivariate stochastic inversion copula
#'
#' The distribution function `C(u1, u2) = P(U1 <= u1, U2 <= u2)` of a
#' \linkS4class{bsicopula}, to go with [dbsicopula()], [hbsicopula()] and
#' [rbsicopula()].
#'
#' With `randomizermod = NULL` the density is `c_V(T1(u1), T2(u2))`
#' ([dbsicopula()]). When both `udp1` and `udp2` are piecewise linear --
#' [vlinear()], [vsymmetric()], [udpzigzag()], [udpcosine()], and shuffles
#' including [udpid()] and [udpflip()] -- the integral is a finite sum. Each
#' linear piece of a transformation is selected, given the carrier value, with
#' a constant probability (the width of the piece for all of these but the
#' shuffles, for which it is `1`), and `C(u1, u2)` is the sum over pairs of
#' pieces of the product of those probabilities and the mass the base copula
#' puts on the rectangle of carrier values the two pieces sweep up to
#' `(u1, u2)`. For two linear v-transforms this is Proposition S3 in the
#' supplementary material of Dias, Han and McNeil. The result is exact to the
#' accuracy of the base copula's CDF.
#'
#' Any other transformation (the non-linear v-transforms [v2p()], [v2b()],
#' [v3p()], [v3b()], and the polynomial and cosine families) is handled by
#' quadrature. The CDF is the integral over one non-linear margin of the
#' h-function of the other, conditioning on the carrier value of the first;
#' each integral is split at the break points of the transformation
#' (where it stops being smooth) and done by tanh-sinh quadrature with
#' `nodes` nodes per cell, which copes with the infinite slopes at cell ends.
#' With one non-linear margin the cost is `nodes` base-copula evaluations per
#' cell and point; with two it is `nodes^2` per pair of cells, so lower
#' `nodes` for large samples. Accuracy is about `1e-8` at the default.
#'
#' With a randomizer model (a \linkS4class{randsdvine} or
#' \linkS4class{randmixture}) the density is `c_V(v1, v2) w(u1, u2)`
#' ([dbsicopula()]), with a weight `w` that depends on which branches of the
#' transformations are selected, so no piece-by-piece closed form exists even
#' for linear transformations. `pbsicopula()` and [hbsicopula()] then
#' integrate [dbsicopula()] itself by the quadrature above, whatever the
#' transformations are, with the CDF a double integral costing `nodes^2`
#' evaluations of the density per pair of cells: use a smaller `nodes`
#' (30 to 40 is usually plenty for a \linkS4class{randsdvine}) for more than
#' a few points. The integrand of a \linkS4class{randmixture} jumps where its
#' selector changes value, which tanh-sinh handles poorly: expect about
#' `1e-3` accuracy there.
#'
#' Exactly on the boundary of the unit square the result is exact; for a
#' `bicop_dist` base copula values elsewhere inherit \pkg{rvinecopulib}'s
#' clipping of its arguments to `[1e-10, 1 - 1e-10]`.
#'
#' @param u1,u2 numeric vectors with values in `[0, 1]`, of equal length or
#'   with a length-1 argument recycled to the length of the other. `u1` may
#'   instead be a two-column matrix with `u2` omitted: `pbsicopula(U, object =
#'   bc)`.
#' @param object an object of class \linkS4class{bsicopula}. Its base copula may be a `bicop_dist` or a
#'   `parCopula` (via [copula::pCopula()]).
#' @param nodes number of tanh-sinh quadrature nodes per cell for
#'   transformations that are not piecewise linear; ignored where the result
#'   is exact.
#'
#' @return A numeric vector of CDF values, one per pair.
#' @references
#' Dias, A., Han, J. and McNeil, A. J. (2026). GARCH copulas, v-transforms and
#' D-vines for stochastic volatility. *Journal of Multivariate Analysis*.
#' @seealso [hbsicopula()], [dbsicopula()], [rbsicopula()].
#' @include bsicopula.R
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bc <- bsicopula(rvinecopulib::bicop_dist("clayton", 90, 2), vlinear(0.4),
#'     udpzigzag(widths = c(3, 4, 3)))
#'   pbsicopula(c(0.2, 0.5, 0.9), c(0.3, 0.5, 0.6), bc)
#'   pbsicopula(0.5, 0.5, bc)
#' }
pbsicopula <- function(u1, u2 = NULL, object, nodes = 101L) {
  a <- bsicopula_args(u1, u2, object, if (missing(object)) NULL else object, range = TRUE)
  check_nodes(nodes)
  if (!is.null(a$object@randomizermod) || is.null(udplinpieces(a$object@udp1)) || is.null(udplinpieces(a$object@udp2))) {
    return(pbsicopula_quad(a$u1, a$u2, a$object, as.integer(nodes)))
  }
  pbsicopula_exact(a$u1, a$u2, a$object)
}

#' h-functions of a bivariate stochastic inversion copula
#'
#' The conditional distribution functions of a \linkS4class{bsicopula} and
#' their inverses, with the conventions of [rvinecopulib::hbicop()]:
#' `cond_var = 1` gives `h(u2 | u1) = P(U2 <= u2 | U1 = u1)`, the partial
#' derivative of the CDF [pbsicopula()] with respect to `u1`, and `cond_var =
#' 2` gives `P(U1 <= u1 | U2 = u2)`. With `inverse = TRUE` the argument that
#' is not conditioned on is a probability level instead, and the function
#' returns the value with that conditional probability: for `cond_var = 1`,
#' `u1` is the conditioning value, `u2` the level and the result the `u2`
#' with `h(u2 | u1)` equal to the level; for `cond_var = 2`, `u2` is the
#' conditioning value, `u1` the level, and the result the `u1`.
#'
#' With `randomizermod = NULL`, conditioning on `U1 = u1` fixes the carrier
#' value `v1 = T1(u1)`, whatever `udp1` is, and the conditional distribution
#' of `U2` is that of the stochastic inverse of `V2 | V1 = v1`. When `udp2`
#' (the margin integrated over) is piecewise linear -- [vlinear()],
#' [vsymmetric()], [udpzigzag()], [udpcosine()] or a shuffle -- `h` is a
#' finite sum of the base copula's h-function at the ends of the carrier
#' intervals the pieces sweep, and its inverse is also closed form: the piece
#' containing the level is found from the cumulative piece probabilities, and
#' one inverse base-copula h-function finishes the job. For two linear
#' v-transforms these are the formulas of Proposition S3 in the supplementary
#' material of Dias, Han and McNeil, which are those used by `tscopula` for
#' vt-D-vines.
#'
#' When `udp2` (the margin integrated over) is not piecewise linear, `h` is
#' the integral of the density `c_V(v1, T2(x))` over `x` up to `u2`, computed
#' cell by cell with tanh-sinh quadrature (`nodes` per cell), and the inverse
#' finds the cell holding the level from the cumulative cell masses and then
#' solves for `u2` by a bracketed Newton iteration, the density being the
#' derivative. This path is several times slower than the exact one, and the
#' inverse slower again.
#'
#' For a `bicop_dist` base copula the base h-functions are those of
#' \pkg{rvinecopulib}, called with the copula's arguments in its own order
#' (never transposed), so non-exchangeable copulas such as the 90 and 270
#' degree rotations are handled correctly. For a `parCopula` they are central
#' differences of [copula::pCopula()] (accurate to about `1e-9` in the
#' interior), and the inverse is found by bisection; \pkg{copula}'s own
#' `cCopula()` is not used, because it conditions only on the first
#' coordinate, has no inverse for rotated copulas, and for a rotation that
#' flips the conditioned coordinate disagrees with the derivative of the CDF
#' (copula 1.1.7). This path is considerably slower than the
#' \pkg{rvinecopulib} one.
#'
#' @inheritParams pbsicopula
#' @param cond_var the conditioning variable, `1` or `2`.
#' @param inverse logical; return the inverse h-function?
#'
#' @return A numeric vector, one value per pair.
#' @references
#' Dias, A., Han, J. and McNeil, A. J. (2026). GARCH copulas, v-transforms and
#' D-vines for stochastic volatility. *Journal of Multivariate Analysis*.
#' @seealso [pbsicopula()], [dbsicopula()], [rbsicopula()].
#' @include bsicopula.R
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 2), vlinear(0.3), vlinear(0.7))
#'   h <- hbsicopula(0.4, c(0.2, 0.6, 0.9), bc)            # P(U2 <= u2 | U1 = 0.4)
#'   h
#'   hbsicopula(0.4, h, bc, inverse = TRUE)                # recovers u2
#'   hbsicopula(c(0.2, 0.6, 0.9), 0.4, bc, cond_var = 2)   # P(U1 <= u1 | U2 = 0.4)
#' }
hbsicopula <- function(u1, u2 = NULL, object, cond_var = 1L, inverse = FALSE, nodes = 101L) {
  a <- bsicopula_args(u1, u2, object, if (missing(object)) NULL else object, range = TRUE)
  check_nodes(nodes)
  if (!(length(cond_var) == 1L && cond_var %in% 1:2)) {
    stop("'cond_var' must be 1 or 2.", call. = FALSE)
  }
  if (!(is.logical(inverse) && length(inverse) == 1L && !is.na(inverse))) {
    stop("'inverse' must be TRUE or FALSE.", call. = FALSE)
  }
  object <- a$object
  # target margin t: the one integrated over (2 when conditioning on U1)
  t <- 3L - as.integer(cond_var)
  u_g <- if (cond_var == 1L) a$u1 else a$u2
  u_t <- if (cond_var == 1L) a$u2 else a$u1
  udp_t <- slot(object, paste0("udp", t))
  Pt <- udplinpieces(udp_t)
  if (is.null(Pt) || !is.null(object@randomizermod)) {
    # the conditioning value is u_g itself when a randomizer makes the density depend on branches
    g <- if (is.null(object@randomizermod)) udptrans(slot(object, paste0("udp", cond_var)), u_g) else u_g
    N <- as.integer(nodes)
    if (inverse) {
      return(h_inverse_quad(u_t, g, udp_t, object, as.integer(cond_var), N))
    }
    return(h_quad(u_t, g, udp_t, object, as.integer(cond_var), N))
  }
  v_g <- udptrans(slot(object, paste0("udp", cond_var)), u_g)
  if (inverse) {
    return(h_inverse_exact(u_t, v_g, Pt, object@basecopula, as.integer(cond_var)))
  }
  v_t <- udptrans(udp_t, u_t)
  h_exact(u_t, v_t, v_g, Pt, object@basecopula, as.integer(cond_var))
}
