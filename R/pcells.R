# Cell probabilities of a bsicopula: the probability that stochastic
# inversion selects, for each margin, the branch of the udp transformation
# lying in a given cell of [0, 1].

# The branches of udp x at each v, located in the cells bounded by 'breaks'.
# udpinverse() gives one column per branch, left-packed: Z is the
# length(v) x (k + 1) matrix of cumulative selection probabilities (0 first,
# 1 from the last branch on), whose consecutive columns bound the randomizer
# interval that selects each branch, and 'cell' the length(v) x k matrix of
# the cell each branch's pre-image lies in. At a turning point the pre-image
# is listed once per branch meeting there; the first of such a coincident
# pair belongs to the cell on the left. A missing branch (trailing NA) gets
# an empty interval, so its cell, set to 1, is immaterial.
cell_branches <- function(x, v, breaks) {
  M <- udpinverse(x, v, prob = TRUE)
  P <- attr(M, "prob")
  present <- !is.na(M)
  P[!present] <- 0
  n <- nrow(M)
  k <- ncol(M)
  Z <- matrix(0, n, k + 1L)
  for (j in seq_len(k)) {
    Z[, j + 1L] <- pmin(Z[, j] + P[, j], 1)
  }
  upper <- Z[, -1L, drop = FALSE]
  upper[!present] <- 1
  upper[cbind(seq_len(n), rowSums(present))] <- 1
  Z[, -1L] <- upper

  nxt <- cbind(M[, -1L, drop = FALSE], NA_real_)
  first_of_pair <- !is.na(nxt) & abs(nxt - M) < 1e-12
  cell <- findInterval(M, breaks, all.inside = TRUE)
  left <- findInterval(M, breaks, left.open = TRUE, all.inside = TRUE)
  cell[first_of_pair] <- left[first_of_pair]
  cell[!present] <- 1L
  dim(cell) <- dim(M)
  list(Z = Z, cell = cell)
}

# Cell probabilities given V1 = v1[i], V2 = v2[i], for every pair i at once:
# a k1 x k2 x length(v1) array. The conditional randomizer CDF on the grid of
# each pair's interval end points, differenced in both directions, gives the
# probability of every pair of branches; those are then added into the cells
# the branches lie in. (The CDF is 0 wherever either argument is 0, so only
# the other grid points are evaluated.)
pcells_pairs <- function(object, v1, v2, breaks1, breaks2) {
  n <- length(v1)
  b1 <- cell_branches(object@udp1, v1, breaks1)
  b2 <- cell_branches(object@udp2, v2, breaks2)
  k1 <- ncol(b1$cell)
  k2 <- ncol(b2$cell)
  Fz <- randomizer_cdf(v1, v2, object@basecopula, object@randomizermod)

  G <- array(0, c(n, k1 + 1L, k2 + 1L))
  G[, -1L, -1L] <- Fz(
    rep(as.vector(b1$Z[, -1L]), times = k2),
    as.vector(b2$Z[, rep(seq_len(k2) + 1L, each = k1)]),
    rep(seq_len(n), times = k1 * k2)
  )
  lo1 <- seq_len(k1)
  lo2 <- seq_len(k2)
  branch_probs <- G[, lo1 + 1L, lo2 + 1L, drop = FALSE] - G[, lo1, lo2 + 1L, drop = FALSE] -
    G[, lo1 + 1L, lo2, drop = FALSE] + G[, lo1, lo2, drop = FALSE]

  out <- array(0, c(length(breaks1) - 1L, length(breaks2) - 1L, n))
  i <- seq_len(n)
  for (a in lo1) {
    for (b in lo2) {
      at <- cbind(b1$cell[, a], b2$cell[, b], i)
      out[at] <- out[at] + branch_probs[, a, b]
    }
  }
  out
}

# Quadrature for integrating over one carrier given the other: for each
# given[i] (the value of V1 if which = 1, of V2 if which = 2), N nodes for
# the other carrier and weights summing to 1, laid out node-fastest.
#   bicop_dist  substitute w = C(other | given), which turns the integral
#               against the conditional density into a plain average over w
#               in (0, 1): nodes are conditional quantiles at the midpoints
#               of N equal cells, weights 1 / N. No singular weight, however
#               unbounded the copula density.
#   parCopula   the copula package has no conditional quantile for an
#               arbitrary conditioning margin, so use fixed nodes, denser
#               towards 0 and 1 (arcsine-spaced), weighted by the copula
#               density there and normalized.
conditional_nodes <- function(basecopula, given, which, N) {
  n <- length(given)
  t <- (seq_len(N) - 0.5) / N
  g <- rep(given, each = N)
  if (is_bicop_dist(basecopula)) {
    w <- rep(t, times = n)
    other <- if (which == 1L) {
      rvinecopulib::hbicop(cbind(g, w), cond_var = 1, family = basecopula, inverse = TRUE)
    } else {
      rvinecopulib::hbicop(cbind(w, g), cond_var = 2, family = basecopula, inverse = TRUE)
    }
    weight <- rep(1 / N, n * N)
  } else {
    other <- rep(stats::qbeta(t, 0.5, 0.5), times = n)
    dens <- if (which == 1L) {
      basecopula_density(g, other, basecopula)
    } else {
      basecopula_density(other, g, basecopula)
    }
    raw <- matrix(dens / stats::dbeta(other, 0.5, 0.5), N, n)
    weight <- as.vector(sweep(raw, 2, colSums(raw), "/"))
  }
  list(given = g, other = other, weight = weight)
}

# Cell probabilities given one carrier only: pcells_pairs() at the quadrature
# nodes of the other carrier, averaged with their weights. A
# k1 x k2 x length(given) array.
pcells_given_one <- function(object, given, which, breaks1, breaks2, N) {
  q <- conditional_nodes(object@basecopula, given, which, N)
  P <- if (which == 1L) {
    pcells_pairs(object, q$given, q$other, breaks1, breaks2)
  } else {
    pcells_pairs(object, q$other, q$given, breaks1, breaks2)
  }
  d <- dim(P)
  W <- matrix(P, d[1] * d[2]) * rep(q$weight, each = d[1] * d[2])
  out <- rowsum(t(W), rep(seq_along(given), each = N), reorder = FALSE)
  array(t(out), c(d[1], d[2], length(given)))
}

#' Cell probabilities of a bivariate stochastic inversion copula
#'
#' Each udp transformation splits `[0, 1]` into cells, and stochastic
#' inversion places `U_i` in one of them: given the carrier value `V_i`, the
#' randomizer `Z_i` selects one pre-image of `V_i`, i.e. one branch of
#' `udp_i`, and `U_i` falls in the cell that branch lies in. `pcells()` gives
#' the joint probability of every pair of cells, `P(U1 in cell l1, U2 in
#' cell l2)`, conditional on both carriers, on one of them, or on neither.
#'
#' **Cells.** With `cells = "monotone"` the cells are the monotone branches
#' of the transformation: the boundaries are the points where it changes
#' direction or jumps -- the fulcrum of a v-transform, the breakpoints of a
#' \linkS4class{udpzigzag} or \linkS4class{udpcosine}, the turning points
#' of the underlying `g` for \linkS4class{udplegendre},
#' \linkS4class{udplegendrebex} and \linkS4class{udpcosinebex}, the strip
#' boundaries of a shuffle. With `cells = "smooth"` they are the intervals
#' on which the transformation is continuously differentiable, which for the
#' three polynomial classes adds the points where another branch starts or
#' stops, and so gives more cells. For every other class the two coincide.
#' Either way a cell holds at most one pre-image of any value. The cell
#' boundaries are returned as the attributes `"breaks1"` and `"breaks2"`.
#'
#' **What is conditioned on.**
#' * `v1` and `v2` both given: the probabilities given `V1 = v1, V2 = v2`.
#'   This is where the randomizer model acts: with `randomizermod = NULL`
#'   the two selections are independent, so the probabilities are the
#'   products of the two marginal selection probabilities; a
#'   \linkS4class{randsdvine} or \linkS4class{randmixture} moves probability
#'   between the cells. To compare, call `pcells()` on the object and on a
#'   copy with `randomizermod` set to `NULL`.
#' * one of them `NULL`: the probabilities given the other alone, the
#'   average of the above over the conditional distribution of the missing
#'   carrier under `basecopula`.
#' * both `NULL`: the unconditional probabilities, which are the mass the
#'   copula itself puts on each rectangle of cells. Their row and column
#'   sums are the cell widths, since `U1` and `U2` are uniform.
#'
#' The averages are computed by quadrature with `ngrid` nodes per
#' integration (so `ngrid^2` points in the unconditional case). For a
#' `bicop_dist` base copula the nodes are conditional quantiles of the
#' missing carrier, which handles unbounded copula densities exactly. For a
#' parCopula base copula (possible when `randomizermod` is `NULL` or a
#' \linkS4class{randmixture}) the \pkg{copula} package offers no conditional
#' quantile for an arbitrary conditioning margin, so fixed nodes weighted by
#' the copula density are used instead. Either way the accuracy is limited
#' by the integrand, which changes abruptly in the missing carrier at the
#' turning values of its udp transformation: expect about three decimal
#' places at the default `ngrid` (more when that transformation is a
#' v-transform, a zigzag or a shuffle), improving roughly in proportion to
#' `ngrid`. The identities above -- row and column sums equal to the cell
#' widths -- give a direct check.
#'
#' At a turning value (a `v` whose pre-image is a turning point) the
#' probabilities are the one-sided limits, as in [udpinverse()]. Exactly at
#' `v = 0` or `v = 1` results for a `bicop_dist` base copula are those at
#' `1e-10` and `1 - 1e-10`, to which \pkg{rvinecopulib} clips its arguments;
#' conditional distributions under the base copula can change quickly that
#' close to the edge.
#'
#' @param object an object of class \linkS4class{bsicopula}.
#' @param v1,v2 numeric vectors with values in `[0, 1]`, or `NULL`; see
#'   Details. When both are given they are paired, recycled to a common
#'   length, unless `grid = TRUE`.
#' @param cells `"monotone"` (the default) or `"smooth"`; see Details.
#' @param grid logical; when both `v1` and `v2` are given, evaluate at every
#'   combination of their values rather than at pairs?
#' @param ngrid number of quadrature nodes per integration when `v1` or `v2`
#'   is `NULL`.
#' @param drop logical; drop the dimensions indexing `v1`/`v2` when they
#'   have extent 1? The two cell dimensions are always kept.
#'
#' @return A numeric array whose first two dimensions index the cells of
#'   `udp1` and of `udp2`. With `v1` and `v2` paired, or one of them `NULL`,
#'   a third dimension indexes the conditioning values; with `grid = TRUE`
#'   the third and fourth index `v1` and `v2`; in the unconditional case
#'   there is a third dimension of extent 1. With `drop = TRUE` such
#'   dimensions of extent 1 are dropped, so a single conditioning value, or
#'   none, gives a matrix. Each slice over the first two dimensions sums to
#'   1. The attributes `"breaks1"` and `"breaks2"` hold the cell boundaries.
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
#' @include bsicopula.R
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   bd <- rvinecopulib::bicop_dist
#'   bc <- bsicopula(bd("gaussian", 0, 0.5), vlinear(0.4), udpzigzag(widths = c(3, 4, 3)),
#'     randsdvine(bd("gaussian", 0, 0.7), bd("clayton", 0, 1.5), bd("gumbel", 0, 1.8))
#'   )
#'   pcells(bc, 0.3, 0.6)
#'
#'   # against independent randomizers
#'   bc0 <- bc
#'   bc0@randomizermod <- NULL
#'   pcells(bc0, 0.3, 0.6)
#'
#'   # given V1 alone, and unconditionally
#'   pcells(bc, v1 = 0.3)
#'   pcells(bc)
#'
#'   # one pair of cells over a grid, for a contour plot
#'   v <- (seq_len(25) - 0.5) / 25
#'   P <- pcells(bc, v, v, grid = TRUE)
#'   contour(v, v, P[2, 2, , ], xlab = "v1", ylab = "v2")
#' }
pcells <- function(object, v1 = NULL, v2 = NULL, cells = c("monotone", "smooth"),
                   grid = FALSE, ngrid = 400L, drop = TRUE) {
  if (!methods::is(object, "bsicopula")) {
    stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
  }
  cells <- match.arg(cells)
  for (arg in c("grid", "drop")) {
    val <- get(arg)
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", arg), call. = FALSE)
    }
  }
  if (!is.numeric(ngrid) || length(ngrid) != 1L || is.na(ngrid) || ngrid < 2 ||
    ngrid != round(ngrid)) {
    stop("'ngrid' must be a single integer of at least 2.", call. = FALSE)
  }
  check_v <- function(v, name) {
    if (is.null(v)) {
      return(NULL)
    }
    if (!is.numeric(v) || !length(v) || anyNA(v) || any(v < 0 | v > 1)) {
      stop(sprintf("'%s' must be NULL or a numeric vector with values in [0, 1].", name),
        call. = FALSE
      )
    }
    as.numeric(v)
  }
  v1 <- check_v(v1, "v1")
  v2 <- check_v(v2, "v2")

  breaksof <- if (cells == "monotone") udpmonobreaks else udpbreaks
  breaks1 <- breaksof(object@udp1)
  breaks2 <- breaksof(object@udp2)

  if (!is.null(v1) && !is.null(v2)) {
    if (grid) {
      out <- pcells_pairs(object, rep(v1, times = length(v2)), rep(v2, each = length(v1)),
        breaks1, breaks2
      )
      dim(out) <- c(dim(out)[1:2], length(v1), length(v2))
    } else {
      n <- max(length(v1), length(v2))
      out <- pcells_pairs(object, rep_len(v1, n), rep_len(v2, n), breaks1, breaks2)
    }
  } else if (!is.null(v1)) {
    out <- pcells_given_one(object, v1, 1L, breaks1, breaks2, ngrid)
  } else if (!is.null(v2)) {
    out <- pcells_given_one(object, v2, 2L, breaks1, breaks2, ngrid)
  } else {
    # unconditional: V1 is uniform, so average the probabilities given V1
    # over the midpoints of ngrid equal cells
    given <- pcells_given_one(object, (seq_len(ngrid) - 0.5) / ngrid, 1L, breaks1, breaks2, ngrid)
    out <- array(rowMeans(given, dims = 2L), c(dim(given)[1:2], 1L))
  }

  d <- dim(out)
  if (drop) {
    d <- c(d[1:2], d[-(1:2)][d[-(1:2)] > 1L])
    dim(out) <- d
  }
  dimnames(out) <- c(
    list(cell1 = as.character(seq_len(d[1])), cell2 = as.character(seq_len(d[2]))),
    rep(list(NULL), length(d) - 2L)
  )
  attr(out, "breaks1") <- breaks1
  attr(out, "breaks2") <- breaks2
  out
}
