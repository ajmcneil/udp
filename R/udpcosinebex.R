# udp transformations built from finite sums of orthonormal cosine functions.
#
# For weights c_1, ..., c_N and the orthonormal cosine functions
# Omega_j(u) = sqrt(2) * (-1)^j * cos(j * pi * u) on [0, 1], let
# g(u) = sum_j c_j Omega_j(u) and let F be the distribution function of g(U)
# for U uniform. The transformation is T(u) = F(g(u)), uniform-distribution-
# preserving for the same reason as udplegendrebex(). udpcosine(j) is the
# one-hot special case coef = c(rep(0, j - 1), 1 / sqrt(2)).
#
# The Chebyshev identity cos(j*theta) = T_j(cos(theta)), with theta = pi*u and
# x = cos(pi*u), turns g into an ordinary polynomial in a DIFFERENT variable:
# h(x) = sqrt(2) * sum_j c_j * (-1)^j * T_j(x), a polynomial on [-1, 1], with
# g(u) = h(cos(pi*u)). Every low-level helper below therefore works in x, not
# u, and results are converted back via u = acos(x) / pi where needed --
# unlike udplegendrebex.R, none of udplegendre.R's [0, 1]-domain helpers
# (legendre_realroots(), legendre_turnpoints()) apply directly, since this
# class's natural domain is [-1, 1]. Only poly_deriv_coef() and polyval()
# (pure coefficient algebra, no domain built in) are reused from there.
#
# Because x = cos(pi*u) is a smooth bijection [0, 1] -> [-1, 1] with
# dx/du != 0 on the OPEN interval (0, 1), turning points of g in u correspond
# exactly to turning points of h in x, and the same panel-boundary logic from
# udplegendrebex.R applies: split panels at h's turning-point critical values
# AND at h(1) = g(0), h(-1) = g(-1)... == g(1) whenever they land strictly
# inside [lbound, ubound]. It matters MORE here than for udplegendrebex: every
# Omega_j has zero derivative at both u = 0 and u = 1 (sin(j*pi*0) =
# sin(j*pi*1) = 0 for every integer j), so g'(0) = g'(1) = 0 identically for
# ANY coefficient choice -- a domain-endpoint crossing is a genuine one-sided
# square-root singularity here, not just an ordinary finite-slope kink.
#
# That same structural fact (g'(0) = g'(1) = 0 always) means udpderiv()'s
# general formula f(g(u)) * g'(u) is wrong at u = 0 and u = 1: it would
# evaluate to 0 there, but the true one-sided derivative is not 0. Working
# from the same local-expansion argument that gives the interior turning-point
# derivative (T's one-sided derivative at a nondegenerate interior turning
# point is exactly +-2 * (coincidence multiplicity), independent of curvature
# magnitude, because two symmetric branches contribute), a boundary point has
# only ONE branch contributing, halving that constant to +-1 per coincidence
# -- and because x'(0) = x'(1) = 0, the controlling second derivative of g at
# the boundary comes entirely from h's FIRST derivative there (the x'' term),
# not h'':
#   g''(0) = h'(1)  * x''(0) = h'(1)  * (-pi^2)  ==>  sign(g''(0)) = -sign(h'(1))
#   g''(1) = h'(-1) * x''(1) = h'(-1) * ( pi^2)  ==>  sign(g''(1)) =  sign(h'(-1))
# For a GENERIC coef this is the whole story (multiplicity 1: T'(0+) =
# -sign(h'(1)), T'(1-) = -sign(h'(-1))), verified against udpcosine(1)'s
# closed form. But a boundary value can ALSO coincide with an interior
# turning point's critical value, or with the other boundary's value -- always
# true of the pure single-term reduction to udpcosine(d), whose equally
# spaced troughs/peaks share every extreme value -- in which case the true
# one-sided derivative is the coincidence count times as large (weight 2 per
# coincident interior turning point, 1 for the other boundary point, matching
# the weighting already used for interior mult); see udpderiv() below.

# Ascending monomial coefficients of the Chebyshev polynomial T_degree(x).
chebyshev_coef <- function(degree) {
  if (degree == 0L) {
    return(1)
  }
  prev <- 1 # T_0
  cur <- c(0, 1) # T_1 = x
  if (degree == 1L) {
    return(cur)
  }
  for (n in seq_len(degree - 1L)) {
    # 2x * cur in ascending coefficients, minus prev (the T_{n+1} recurrence)
    shifted <- c(0, 2 * cur)
    prev <- c(prev, numeric(length(shifted) - length(prev)))
    nxt <- shifted - prev
    prev <- cur
    cur <- nxt
  }
  cur
}

# Ascending monomial coefficients of h(x) = sum_j coef[j] * sqrt(2) * (-1)^j *
# T_j(x), the polynomial g(u) = sum_j coef[j] * Omega_j(u) becomes under
# x = cos(pi * u).
chebyshev_bex_coef <- function(coef) {
  n <- length(coef)
  cfs <- numeric(n + 1L)
  for (j in seq_len(n)) {
    if (coef[j] == 0) next
    cj <- chebyshev_coef(j) * (sqrt(2) * (-1)^j * coef[j])
    cfs[seq_along(cj)] <- cfs[seq_along(cj)] + cj
  }
  cfs
}

# Real roots in [-1, 1] of h(x) = y. Identical in method to
# legendre_realroots() (udplegendre.R) but searches [-1, 1] instead of [0, 1].
chebyshev_realroots <- function(coef, y, tol = 1e-7) {
  cc <- coef
  cc[1] <- cc[1] - y
  z <- polyroot(cc)
  re <- Re(z)
  cand <- re[abs(Im(z)) < 1e-6 * pmax(1, abs(re)) & re > -1 - tol & re < 1 + tol]
  if (!length(cand)) {
    return(numeric(0))
  }
  cand <- pmin(pmax(cand, -1), 1)
  resid <- abs(polyval(coef, cand) - y)
  slope <- pmax(abs(polyval(poly_deriv_coef(coef), cand)), 1)
  cand <- sort(cand[resid < 1e-6 * slope + 1e-9])
  if (length(cand) > 1L) {
    cand <- cand[c(TRUE, diff(cand) > tol)]
  }
  cand
}

# Interior turning points of h (roots of h' in (-1, 1)). Identical in method
# to legendre_turnpoints() but searches (-1, 1) instead of (0, 1).
chebyshev_turnpoints <- function(coefD) {
  if (length(coefD) < 2L) {
    return(numeric(0))
  }
  z <- polyroot(coefD)
  tp <- Re(z)[abs(Im(z)) < 1e-6]
  sort(tp[tp > -1 & tp < 1])
}

# Range of h on [-1, 1]: the minimum and maximum over the candidate set of the
# two endpoints and the interior turning points.
chebyshev_bounds <- function(coef, coefD) {
  cand <- c(-1, 1, chebyshev_turnpoints(coefD))
  vals <- polyval(coef, cand)
  c(min(vals), max(vals))
}

# Exact F(y) = |{u in [0, 1] : g(u) <= y}| for a single y, where
# g(u) = h(cos(pi * u)). Unlike legendre_measure_bounded() (udplegendrebex.R),
# this cannot count Lebesgue measure directly in x: x = cos(pi * U) is NOT
# uniform when U is (it has the arcsine distribution), so the roots of
# h(x) = y are found in x, converted to u via u = acos(x) / pi, and the
# breaks-and-midpoint measure is taken in u, which IS the right space (U is
# uniform).
chebyshev_measure_bounded <- function(coef, y, lbound, ubound) {
  if (y <= lbound) {
    return(0)
  }
  if (y >= ubound) {
    return(1)
  }
  rx <- chebyshev_realroots(coef, y)
  ru <- sort(acos(pmin(pmax(rx, -1), 1)) / pi)
  breaks <- c(0, ru, 1)
  mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  gmids <- polyval(coef, cos(pi * mids))
  sum(diff(breaks)[gmids <= y])
}

# g'(u) = h'(x) * dx/du, x = cos(pi * u), dx/du = -pi * sin(pi * u). The chain
# rule for the derivative of a udpcosinebex transformation's underlying g;
# used by udpinverse()'s pre-image weights and udpderiv()'s interior formula.
cosine_sum_gderiv <- function(cfsD, u) {
  -pi * sin(pi * u) * polyval(cfsD, cos(pi * u))
}

#' Class of finite orthonormal-cosine-expansion udp transformations
#'
#' A udpcosinebex transformation is `T(u) = F(g(u))`, the composition of
#' `g(u) = sum_j coef[j] * Omega_j(u)` -- a finite expansion in the orthonormal
#' cosine functions `Omega_j(u) = sqrt(2) * (-1)^j * cos(j * pi * u)` on
#' `[0, 1]` -- with the distribution function `F` of `g(U)` for `U` uniform.
#' It is a uniform-distribution-preserving map of `[0, 1]`: for `U` uniform,
#' `udptrans(x, U)` is again uniform. \linkS4class{udpcosine} is the one-term
#' special case; see [udpcosinebex()] for how the two relate.
#'
#' The identity `cos(j*theta) = T_j(cos(theta))` turns `g` into an ordinary
#' polynomial `h(x) = sqrt(2) * sum_j coef[j] * (-1)^j * T_j(x)` in the
#' DIFFERENT variable `x = cos(pi * u)`, with `g(u) = h(cos(pi * u))`. `F` and
#' its inverse are built exactly as for \linkS4class{udplegendrebex} --
#' exact real roots from polynomial root finding, panel splines of `F` and
#' `F^{-1}` between consecutive turning-point critical values, regularized
#' against the square-root branch point at each -- but working in `x`, with
#' results converted back to `u` via `u = acos(x) / pi`.
#'
#' @slot coef numeric; the weights on `Omega_1, ..., Omega_degree`, trailing
#'   zeros trimmed.
#' @slot degree integer; the index of the last nonzero weight.
#' @slot cfs,cfsD ascending monomial coefficients of `h` (a polynomial in `x`,
#'   not `u`) and of its derivative.
#' @slot lbound,ubound the range of `g` (equivalently of `h`) on `[0, 1]`
#'   (equivalently on `[-1, 1]`).
#' @slot Tfun,Qfun functions evaluating `T` and `F^{-1}`.
#'
#' @seealso [udpcosinebex()] to construct one, [udptrans()] to evaluate it.
#' @include udp-package.R udplegendre.R
#' @export
#'
#' @references
#' McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
#' of non-monotonic dependence. \href{https://arxiv.org/abs/2512.10828}{arXiv:2512.10828}
setClass("udpcosinebex",
  contains = "udp",
  slots = list(
    coef = "numeric", degree = "integer", cfs = "numeric", cfsD = "numeric",
    lbound = "numeric", ubound = "numeric", Tfun = "function", Qfun = "function"
  )
)

#' Construct a finite orthonormal-cosine-expansion udp transformation
#'
#' @param coef a numeric vector of weights, `coef[j]` multiplying the
#'   orthonormal cosine function `Omega_j(u) = sqrt(2) * (-1)^j * cos(j * pi *
#'   u)`. Trailing zeros are trimmed to find the effective degree; at least
#'   one entry must be nonzero. There is no `j = 0` (constant) term, for the
#'   same reason as \linkS4class{udplegendrebex}.
#' @param ngrid number of grid points (at least 3) at which `F` is evaluated
#'   exactly before interpolation.
#'
#' @return An object of class \linkS4class{udpcosinebex}.
#'
#' @details
#' A one-hot `coef` reproduces the corresponding \linkS4class{udpcosine}
#' exactly, up to the `1 / sqrt(2)` rescaling needed to undo the orthonormal
#' `sqrt(2)` factor: `udptrans(udpcosinebex(c(0, 0, 1 / sqrt(2))), u)` equals
#' `udptrans(udpcosine(3), u)`. As for \linkS4class{udplegendrebex}, flipping
#' the sign of every entry of `coef` replaces `T` with its complement, `1 - T`.
#'
#' @export
#'
#' @examples
#' udpcosinebex(c(0.44, -0.35, 0.56, -0.21, 0.36))
#' plot(udpcosinebex(c(0.44, -0.35, 0.56, -0.21, 0.36)), embellish = "colour")
udpcosinebex <- function(coef, ngrid = 513L) {
  if (!is.numeric(coef) || length(coef) < 1L || anyNA(coef)) {
    stop("'coef' must be a numeric vector with no missing values.", call. = FALSE)
  }
  nz <- which(coef != 0)
  if (!length(nz)) {
    stop("'coef' must have at least one nonzero entry.", call. = FALSE)
  }
  degree <- max(nz)
  coef <- as.numeric(coef[seq_len(degree)])
  if (degree > 12L) {
    warning(
      "udpcosinebex(): effective degree > 12 is numerically unreliable ",
      "(polyroot on large alternating coefficients).",
      call. = FALSE
    )
  }
  if (!is.numeric(ngrid) || length(ngrid) != 1L || is.na(ngrid) || ngrid < 3) {
    stop("'ngrid' must be a single number of at least 3.", call. = FALSE)
  }
  ngrid <- as.integer(ngrid)

  cfs <- chebyshev_bex_coef(coef)
  cfsD <- poly_deriv_coef(cfs)
  bounds <- chebyshev_bounds(cfs, cfsD)
  lbound <- bounds[1L]
  ubound <- bounds[2L]

  # Panel boundaries: turning-point critical values, plus h(1) = g(0) and
  # h(-1) = g(1) whenever either lies strictly inside the range -- see the
  # file header for why this matters more here than for udplegendrebex.
  yv <- sort(c(polyval(cfs, c(-1, 1)), polyval(cfs, chebyshev_turnpoints(cfsD))))
  if (length(yv)) {
    yv <- yv[c(TRUE, diff(yv) > 1e-7)]
  }
  ye <- yv[yv > lbound + 1e-7 & yv < ubound - 1e-7]
  ypanels <- c(lbound, ye, ubound)
  npan <- length(ypanels) - 1L
  per <- max(40L, ceiling(ngrid / npan))

  panel <- lapply(seq_len(npan), function(p) {
    a <- ypanels[p]
    b <- ypanels[p + 1L]
    mid <- (a + b) / 2
    half <- (b - a) / 2
    th <- seq(0, pi, length.out = per)
    yy <- mid - half * cos(th)
    yy[c(1L, per)] <- c(a, b)
    fv <- cummax(vapply(yy, chebyshev_measure_bounded, numeric(1),
      coef = cfs, lbound = lbound, ubound = ubound
    ))
    kp <- c(TRUE, diff(fv) > 0)
    list(
      mid = mid, half = half, vhi = fv[per],
      Fofth = stats::splinefun(th, fv, method = "monoH.FC"),
      thofF = stats::splinefun(fv[kp], th[kp], method = "monoH.FC")
    )
  })
  vpanels <- c(0, vapply(panel, `[[`, numeric(1), "vhi"))
  vpanels[npan + 1L] <- 1

  panel_apply <- function(z, breaks, f) {
    p <- pmin(pmax(findInterval(z, breaks, rightmost.closed = TRUE), 1L), npan)
    out <- numeric(length(z))
    for (pp in seq_len(npan)) {
      take <- p == pp
      if (any(take)) out[take] <- f(panel[[pp]], z[take])
    }
    out
  }
  Ffun <- function(y) {
    pmin(pmax(panel_apply(pmin(pmax(y, lbound), ubound), ypanels, function(pl, yy) {
      pl$Fofth(acos(pmin(pmax((pl$mid - yy) / pl$half, -1), 1)))
    }), 0), 1)
  }
  Qfun <- function(v) {
    pmin(pmax(panel_apply(pmin(pmax(v, 0), 1), vpanels, function(pl, vv) {
      pl$mid - pl$half * cos(pmin(pmax(pl$thofF(vv), 0), pi))
    }), lbound), ubound)
  }
  Tfun <- function(u) Ffun(polyval(cfs, cos(pi * pmin(pmax(u, 0), 1))))

  new("udpcosinebex",
    coef = coef, degree = as.integer(degree), cfs = cfs, cfsD = cfsD,
    lbound = lbound, ubound = ubound, Tfun = Tfun, Qfun = Qfun
  )
}

#' @describeIn udptrans Evaluate a udpcosinebex transformation.
#' @export
setMethod("udptrans", "udpcosinebex", function(x, u) {
  out <- pmin(pmax(x@Tfun(pmin(pmax(as.numeric(u), 0), 1)), 0), 1)
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' @describeIn udpinverse Pre-images of a udpcosinebex transformation: a
#'   matrix with `degree` columns holding, for each `v`, the roots in
#'   `[0, 1]` of `g(u) = F^{-1}(v)`, sorted ascending and left-packed with
#'   trailing `NA`. With `prob = TRUE` the `"prob"` attribute weights each
#'   root by `1 / |g'(u)|`, normalized over the row.
#' @export
setMethod("udpinverse", "udpcosinebex", function(x, v, prob = FALSE, ...) {
  vv <- as.numeric(v)
  if (anyNA(vv) || any(vv < 0 | vv > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  k <- x@degree
  y <- pmin(pmax(x@Qfun(vv), x@lbound), x@ubound)
  roots <- lapply(y, function(yi) {
    rx <- chebyshev_realroots(x@cfs, yi)
    ru <- sort(acos(pmin(pmax(rx, -1), 1)) / pi)
    if (length(ru) > k) ru[seq_len(k)] else ru
  })
  lens <- lengths(roots)
  M <- matrix(NA_real_, length(vv), k)
  M[cbind(rep(seq_along(roots), lens), sequence(lens))] <-
    unlist(roots, use.names = FALSE)
  if (prob) {
    present <- !is.na(M)
    w <- matrix(0, nrow(M), k)
    w[present] <- 1 / abs(cosine_sum_gderiv(x@cfsD, M[present]))
    attr(M, "prob") <- finalise_prob(w, present)
  }
  M
})

#' @describeIn udpderiv `T'(u) = f(g(u)) * g'(u)` away from turning points and
#'   the boundary, `f` the density of `g(U)` as in \linkS4class{udplegendrebex}.
#'   At an interior turning point of `g` this is `2 * m` or `-2 * m` according
#'   to the sign of `g''`, exactly as for \linkS4class{udplegendrebex}. At
#'   `u = 0` or `u = 1` the general formula does not apply: `g'` is
#'   identically `0` at both (`sin(j * pi * 0) = sin(j * pi * 1) = 0` for
#'   every `j`, so this holds for any `coef`), but `T` still has a one-sided
#'   derivative there, `m0 * -sign(h'(1))` at `u = 0` and `m1 * -sign(h'(-1))`
#'   at `u = 1`, where `m0`/`m1` count coincidences of `g(0)`/`g(1)` with an
#'   interior turning-point value (weight `2`) or with each other (weight
#'   `1`) -- see the file header for the derivation.
#' @export
setMethod("udpderiv", "udpcosinebex", function(x, u) {
  uu <- pmin(pmax(as.numeric(u), 0), 1)
  cfs <- x@cfs
  cfsD <- x@cfsD
  tp_x <- chebyshev_turnpoints(cfsD)
  cfsDD <- if (length(tp_x)) poly_deriv_coef(cfsD) else numeric(0)
  tp_u <- if (length(tp_x)) sort(acos(pmin(pmax(tp_x, -1), 1)) / pi) else numeric(0)
  tpval <- if (length(tp_x)) polyval(cfs, tp_x) else numeric(0)
  mult <- vapply(tpval, function(v) sum(abs(tpval - v) < 1e-6), integer(1))
  tol <- 1e-5
  hp1 <- polyval(cfsD, 1)
  hpm1 <- polyval(cfsD, -1)
  # Boundary multiplicity: g(0) = h(1) and g(1) = h(-1) can coincide with an
  # interior turning point's critical value (weight 2, as for any coincident
  # turning point) or with each other (weight 1) -- always true, for
  # instance, of the pure single-term reduction to udpcosine(d), whose
  # equally-spaced troughs/peaks share every extreme value. Missing this
  # multiplies the true one-sided derivative by 1/mult.
  y0 <- polyval(cfs, 1)
  y1 <- polyval(cfs, -1)
  mult0 <- 1L + 2L * sum(abs(tpval - y0) < 1e-6) + as.integer(abs(y1 - y0) < 1e-6)
  mult1 <- 1L + 2L * sum(abs(tpval - y1) < 1e-6) + as.integer(abs(y0 - y1) < 1e-6)
  out <- vapply(uu, function(ui) {
    if (ui < tol) {
      return(mult0 * (-sign(hp1)))
    }
    if (ui > 1 - tol) {
      return(mult1 * (-sign(hpm1)))
    }
    if (length(tp_u)) {
      j <- which.min(abs(ui - tp_u))
      if (abs(ui - tp_u[j]) < tol) {
        return(-2 * mult[j] * sign(polyval(cfsDD, tp_x[j])))
      }
    }
    y <- polyval(cfs, cos(pi * ui))
    rx <- chebyshev_realroots(cfs, y)
    ru <- acos(pmin(pmax(rx, -1), 1)) / pi
    fY <- sum(1 / abs(cosine_sum_gderiv(cfsD, ru)))
    fY * cosine_sum_gderiv(cfsD, ui)
  }, numeric(1))
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

# Breakpoints: 0, 1, and every non-smooth point of T. Similar in method to
# udplegendrebex's, but with the same broadened critical-value set (see
# there): also search for other pre-images of g(0) = h(1) and g(1) = h(-1),
# not just of the turning-point values. u = 0 and u = 1 are always included
# regardless (they always are, unconditionally, for this class); if tp_x is
# empty, g is injective (no turning point means no repeated value), so no
# cross points of either kind are possible and the early return is exact.
setMethod("udpbreaks", "udpcosinebex", function(x) {
  cfs <- x@cfs
  cfsD <- x@cfsD
  tp_x <- chebyshev_turnpoints(cfsD)
  if (!length(tp_x)) {
    return(c(0, 1))
  }
  tol <- 5e-5
  tp_u <- sort(acos(pmin(pmax(tp_x, -1), 1)) / pi)
  yv <- sort(c(polyval(cfs, tp_x), polyval(cfs, c(-1, 1))))
  yv <- yv[c(TRUE, diff(yv) > 1e-7)]
  base <- c(0, 1, tp_u)
  cross <- unlist(lapply(yv, function(y) {
    rx <- chebyshev_realroots(cfs, y)
    ru <- acos(pmin(pmax(rx, -1), 1)) / pi
    ru[vapply(ru, function(ri) all(abs(ri - base) > tol), logical(1))]
  }))
  pts <- sort(c(base, cross))
  pts[c(TRUE, diff(pts) > tol)]
})

#' @describeIn pcoincide Integrate `sum_j p_j(v)^2` as in the default method,
#'   but split the range at the images of the turning points of `g`, where the
#'   integrand has a corner, so each piece is smooth.
#' @export
setMethod("pcoincide", "udpcosinebex", function(x) {
  tp_x <- chebyshev_turnpoints(x@cfsD)
  tp_u <- if (length(tp_x)) sort(acos(pmin(pmax(tp_x, -1), 1)) / pi) else numeric(0)
  breaks <- pmin(pmax(x@Tfun(tp_u), 0), 1)
  integrate_collision(x, breaks)
})

#' Plot method for the udpcosinebex class
#'
#' Draws the graph of the udpcosinebex transformation over thin gridlines
#' marking its A-partition and T-partition: vertical lines at the points of
#' `udpbreaks()` (the A-partition, `0` and `1` excluded as redundant with the
#' plot's own border), the coarsest partition of `[0, 1]` into intervals on
#' which `T` is continuously differentiable; horizontal lines at the distinct
#' values `T` takes at those points (the T-partition), `0` and `1` again
#' excluded. A point with a one-sided vertical tangent but no fold reaches
#' the same `T` value as the turning point whose critical value it shares, so
#' it contributes a vertical line without adding a new horizontal one; `u = 0`
#' and `u = 1` are always A-partition members but need not be trivial
#' T-partition values (every `Omega_j` has zero derivative at both, so they
#' behave like one-sided turning points -- see the file header), so they are
#' handled separately in each role.
#'
#' @param x an object of class \linkS4class{udpcosinebex}.
#' @param n number of points at which to evaluate the transformation.
#' @param xlab,ylab axis labels.
#' @param embellish style of the gridlines: `"none"` (the default) to omit
#'   them, `"colour"` for red, or `"bw"` for grey.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(udpcosinebex(c(0.44, -0.35, 0.56, -0.21, 0.36)), embellish = "colour")
#' plot(udpcosinebex(c(0, 0.6, 0, -0.4)), embellish = "bw")
setMethod("plot", c(x = "udpcosinebex", y = "missing"),
  function(x, n = 500L, xlab = "u", ylab = "T(u)",
           embellish = c("none", "colour", "bw"), ...) {
    emb <- plot_embellish(embellish)
    b <- udpbreaks(x)
    u <- sort(unique(c(seq(0, 1, length.out = n), b)))
    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", asp = 1,
      xlab = xlab, ylab = ylab, ...
    )
    if (!is.null(emb)) {
      # A-partition: dropping 0/1 only omits a line redundant with the plot's
      # own border, not a T-partition value -- T(0)/T(1) are computed below
      # from the full, unfiltered `b`.
      bu <- b[b > 0 & b < 1]
      if (length(bu)) segments(bu, 0, bu, 1, col = emb$grid, lwd = 0.5)
      # T-partition: T at every A-partition point, including the domain
      # edges (always members by definition, and not necessarily trivial --
      # see udpbreaks()); values at/near 0 or 1 are dropped as redundant
      # with the border, not because their point was.
      # Two break points sharing a critical value give T() the same value
      # exactly in principle, but each is evaluated through the panel spline
      # independently, so the results can differ at the interpolation-error
      # level (down to ~1e-7 seen in practice, far below any genuine gap
      # between distinct T-partition values) -- deduplicate by tolerance, not
      # exact equality, matching the pattern used for y-values elsewhere in
      # this file.
      tv <- sort(udptrans(x, b))
      tv <- tv[c(TRUE, diff(tv) > 1e-5)]
      tv <- tv[tv > 1e-9 & tv < 1 - 1e-9]
      if (length(tv)) segments(0, tv, 1, tv, col = emb$grid, lwd = 0.5)
    }
    lines(u, udptrans(x, u), lwd = 1.5)
  }
)
