# udp transformations built from finite sums of orthonormal shifted Legendre
# polynomials.
#
# For weights c_1, ..., c_N and the orthonormal shifted Legendre polynomials
# P_j(u) = sqrt(2j + 1) L_j(u) on [0, 1], let g(u) = sum_j c_j P_j(u) and let F
# be the distribution function of g(U) for U uniform. The transformation is
# T(u) = F(g(u)), uniform-distribution-preserving for the same reason as
# udplegendre(): F is strictly increasing on the range of g, so T has exactly
# the (up to N) pre-images of g, and the stochastic-inverse weights 1 / |T'|
# reduce to 1 / |g'| after normalization (the common F' factor cancels).
# udplegendre(j) is the one-hot special case coef = c(rep(0, j - 1), 1), up to
# the sign of coef[j] -- see the constructor documentation.
#
# Unlike a single shifted Legendre polynomial, a general combination g need
# not attain any particular value at u = 0 or u = 1, so (unlike udplegendre.R)
# the range [lbound, ubound] of g has two ends that both have to be found, and
# every place udplegendre.R could hard-code the upper end of the range as 1
# needs the actual ubound instead. legendre_bounds() and
# legendre_measure_bounded() below are what changes; everything else --
# slegendre_coef(), poly_deriv_coef(), polyval(), legendre_realroots(),
# legendre_turnpoints() -- is already written generically over "a polynomial's
# ascending coefficients" and is reused from udplegendre.R unmodified.

# Ascending monomial coefficients of g(u) = sum_j coef[j] * P_j(u), ready to be
# combined with the low-level polynomial helpers in udplegendre.R.
legendre_bex_coef <- function(coef) {
  n <- length(coef)
  cfs <- numeric(n + 1L)
  for (j in seq_len(n)) {
    if (coef[j] == 0) next
    cj <- slegendre_coef(j) * (sqrt(2 * j + 1) * coef[j])
    cfs[seq_along(cj)] <- cfs[seq_along(cj)] + cj
  }
  cfs
}

# Range of g on [0, 1]: the minimum and maximum of g over the candidate set of
# the two endpoints and the interior turning points, since (unlike a pure
# L_j) a general combination need not attain its extremes at u = 0 or u = 1.
legendre_bounds <- function(coef, coefD) {
  cand <- c(0, 1, legendre_turnpoints(coefD))
  vals <- polyval(coef, cand)
  c(min(vals), max(vals))
}

# Exact F(y) = |{u in [0, 1] : g(u) <= y}| for a single y, where g ranges over
# [lbound, ubound]. Identical in method to legendre_measure() (udplegendre.R)
# but tested against the actual ubound of g rather than the literal 1 -- only
# a pure single-term Legendre polynomial is guaranteed to attain 1.
legendre_measure_bounded <- function(coef, y, lbound, ubound) {
  if (y <= lbound) {
    return(0)
  }
  if (y >= ubound) {
    return(1)
  }
  breaks <- c(0, legendre_realroots(coef, y), 1)
  mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  sum(diff(breaks)[polyval(coef, mids) <= y])
}

#' Class of finite orthonormal-Legendre-expansion udp transformations
#'
#' A udplegendrebex transformation is `T(u) = F(g(u))`, the composition of
#' `g(u) = sum_j coef[j] * P_j(u)` -- a finite expansion in the orthonormal
#' shifted Legendre polynomials `P_j(u) = sqrt(2j + 1) L_j(u)` on `[0, 1]` --
#' with the distribution function `F` of `g(U)` for `U` uniform. It is a
#' uniform-distribution-preserving map of `[0, 1]`: for `U` uniform,
#' `udptrans(x, U)` is again uniform. \linkS4class{udplegendre} is the
#' one-term special case; see [udplegendrebex()] for how the two relate.
#'
#' `F` and its inverse are constructed exactly as for \linkS4class{udplegendre}
#' -- exact real roots of `g(u) = y` from [udpinverse()]'s polynomial root
#' finder, panel splines of `F` and `F^{-1}` between consecutive turning-point
#' values of `g`, regularized against the square-root branch point at each --
#' except that the range of `g` need not be `[lbound, 1]`, so both ends of the
#' range are carried as separate slots.
#'
#' @slot coef numeric; the weights on `P_1, ..., P_degree`, trailing zeros
#'   trimmed.
#' @slot degree integer; the index of the last nonzero weight.
#' @slot cfs,cfsD ascending monomial coefficients of `g` and of its
#'   derivative.
#' @slot lbound,ubound the range of `g` on `[0, 1]`.
#' @slot Tfun,Qfun functions evaluating `T` and `F^{-1}`.
#'
#' @seealso [udplegendrebex()] to construct one, [udptrans()] to evaluate it.
#' @include udp-package.R udplegendre.R
#' @export
#'
#' @references
#' McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
#' of non-monotonic dependence. \href{https://arxiv.org/abs/2512.10828}{arXiv:2512.10828}
setClass("udplegendrebex",
  contains = "udp",
  slots = list(
    coef = "numeric", degree = "integer", cfs = "numeric", cfsD = "numeric",
    lbound = "numeric", ubound = "numeric", Tfun = "function", Qfun = "function"
  )
)

#' Construct a finite orthonormal-Legendre-expansion udp transformation
#'
#' @param coef a numeric vector of weights, `coef[j]` multiplying the
#'   orthonormal shifted Legendre polynomial `P_j(u) = sqrt(2j + 1) L_j(u)`.
#'   Trailing zeros are trimmed to find the effective degree; at least one
#'   entry must be nonzero. There is no `j = 0` (constant) term: `T(u) =
#'   F(g(u))` is unchanged by an additive shift of `g`, so a constant term
#'   would carry no information.
#' @param ngrid number of grid points (at least 3) at which `F` is evaluated
#'   exactly before interpolation.
#'
#' @return An object of class \linkS4class{udplegendrebex}.
#'
#' @details
#' Because `P_j` is a positive rescaling of `L_j`, and `T(u) = F_{g(U)}(g(u))`
#' is unchanged by rescaling `g` by any positive constant, a one-hot `coef`
#' with a single positive entry reproduces the corresponding
#' \linkS4class{udplegendre} exactly: `udptrans(udplegendrebex(c(0, 0, 1)),
#' u)` equals `udptrans(udplegendre(3), u)`. Flipping the sign of every entry
#' of `coef` (equivalently, negating `g`) replaces `T` with its complement,
#' `1 - T`.
#'
#' Unlike \linkS4class{udplegendre}, no closed form is used for any effective
#' degree, including 1 or 2: every `coef` goes through the same panel-spline
#' construction, for [udplegendre()] itself if the exact, fast single-term
#' closed forms are wanted.
#'
#' Panel boundaries are placed at the critical values of `g`'s turning points,
#' as for \linkS4class{udplegendre}, and additionally at `g(0)` and `g(1)`
#' whenever either lies strictly inside the range: a pure `L_j` always
#' attains `lbound` or `1` exactly at `u = 0` or `u = 1`, but a general
#' combination need not, and wherever `g(0)` or `g(1)` sits strictly between
#' `lbound` and `ubound`, `F` has a genuine (finite-slope) kink there -- a
#' second branch starts or stops contributing at that value. Splitting the
#' panel there, rather than asking one spline to interpolate through the
#' kink, is what keeps accuracy close to \linkS4class{udplegendre}'s.
#'
#' @export
#'
#' @examples
#' udplegendrebex(c(0.44, -0.35, 0.56, -0.21, 0.36))
#' plot(udplegendrebex(c(0.44, -0.35, 0.56, -0.21, 0.36)), embellish = "colour")
udplegendrebex <- function(coef, ngrid = 513L) {
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
      "udplegendrebex(): effective degree > 12 is numerically unreliable ",
      "(polyroot on large alternating coefficients).",
      call. = FALSE
    )
  }
  if (!is.numeric(ngrid) || length(ngrid) != 1L || is.na(ngrid) || ngrid < 3) {
    stop("'ngrid' must be a single number of at least 3.", call. = FALSE)
  }
  ngrid <- as.integer(ngrid)

  cfs <- legendre_bex_coef(coef)
  cfsD <- poly_deriv_coef(cfs)
  bounds <- legendre_bounds(cfs, cfsD)
  lbound <- bounds[1L]
  ubound <- bounds[2L]

  # Interior panel-boundary values: the critical values of g's turning points
  # (where F has a square-root singularity, as in udplegendre.R) AND g(0),
  # g(1) (where F has an ordinary finite-slope kink instead -- unlike a pure
  # L_j, whose value at u = 0 and u = 1 always equals lbound or ubound
  # exactly, a general combination can have g(0) or g(1) strictly inside the
  # range, meaning a second branch starts or stops contributing to F at that
  # value without either endpoint being flagged as a break. Left unsplit, the
  # panel spline is asked to interpolate straight through a real kink; adding
  # it as a panel boundary is what lets each side be smooth. (Unlike a pure
  # L_j, a combination can also be globally monotonic at any degree, so
  # turning points -- and hence yv -- can be empty; indexing an empty vector
  # with a length-1 logical, as udplegendre.R's version of this line does,
  # returns NA rather than numeric(0), so the empty case is guarded.)
  yv <- sort(c(polyval(cfs, c(0, 1)), polyval(cfs, legendre_turnpoints(cfsD))))
  if (length(yv)) {
    yv <- yv[c(TRUE, diff(yv) > 1e-7)]
  }
  ye <- yv[yv > lbound + 1e-7 & yv < ubound - 1e-7]
  ypanels <- c(lbound, ye, ubound)
  npan <- length(ypanels) - 1L
  per <- max(40L, ceiling(ngrid / npan))

  # Same square-root-branch-point substitution as udplegendre(): on panel
  # p = [a, b] between consecutive extreme values, y = mid - half * cos(theta)
  # regularizes F near both ends of the panel.
  panel <- lapply(seq_len(npan), function(p) {
    a <- ypanels[p]
    b <- ypanels[p + 1L]
    mid <- (a + b) / 2
    half <- (b - a) / 2
    th <- seq(0, pi, length.out = per)
    yy <- mid - half * cos(th)
    yy[c(1L, per)] <- c(a, b)
    fv <- cummax(vapply(yy, legendre_measure_bounded, numeric(1),
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
  Tfun <- function(u) Ffun(polyval(cfs, pmin(pmax(u, 0), 1)))

  new("udplegendrebex",
    coef = coef, degree = as.integer(degree), cfs = cfs, cfsD = cfsD,
    lbound = lbound, ubound = ubound, Tfun = Tfun, Qfun = Qfun
  )
}

#' @describeIn udptrans Evaluate a udplegendrebex transformation.
#' @export
setMethod("udptrans", "udplegendrebex", function(x, u) {
  out <- pmin(pmax(x@Tfun(pmin(pmax(as.numeric(u), 0), 1)), 0), 1)
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' @describeIn udpinverse Pre-images of a udplegendrebex transformation: a
#'   matrix with `degree` columns holding, for each `v`, the roots in
#'   `[0, 1]` of `g(u) = F^{-1}(v)`, sorted ascending and left-packed with
#'   trailing `NA`. With `prob = TRUE` the `"prob"` attribute weights each
#'   root by `1 / |g'(u)|`, normalized over the row.
#' @export
setMethod("udpinverse", "udplegendrebex", function(x, v, prob = FALSE, ...) {
  vv <- as.numeric(v)
  if (anyNA(vv) || any(vv < 0 | vv > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  k <- x@degree
  y <- pmin(pmax(x@Qfun(vv), x@lbound), x@ubound)
  roots <- lapply(y, function(yi) {
    r <- legendre_realroots(x@cfs, yi)
    if (length(r) > k) r[seq_len(k)] else r
  })
  lens <- lengths(roots)
  M <- matrix(NA_real_, length(vv), k)
  M[cbind(rep(seq_along(roots), lens), sequence(lens))] <-
    unlist(roots, use.names = FALSE)
  if (prob) {
    present <- !is.na(M)
    w <- matrix(0, nrow(M), k)
    w[present] <- 1 / abs(polyval(x@cfsD, M[present]))
    attr(M, "prob") <- finalise_prob(w, present)
  }
  M
})

#' @describeIn udpderiv `T'(u) = f(g(u)) * g'(u)`, where `f` is the density of
#'   `g(U)`, `sum(1 / abs(g'(u_i)))` over the pre-images of `g(u)` found by the
#'   same root-finding as [udpinverse()]. At a turning point of `g` this is
#'   `0 * Inf`; it is replaced there by the exact left derivative, `2 * m` or
#'   `-2 * m` according to the sign of `g''`, where `m` is the number of
#'   turning points sharing that critical value (a coincidence, not a
#'   structural feature, for a generic combination). At any other pre-image of
#'   a turning-point value `T` has a one-sided vertical tangent, which the
#'   formula already returns as `Inf` or `-Inf` without a special case.
#' @export
setMethod("udpderiv", "udplegendrebex", function(x, u) {
  uu <- pmin(pmax(as.numeric(u), 0), 1)
  cfs <- x@cfs
  cfsD <- x@cfsD
  tp <- legendre_turnpoints(cfsD)
  cfsDD <- if (length(tp)) poly_deriv_coef(cfsD) else numeric(0)
  tpval <- if (length(tp)) polyval(cfs, tp) else numeric(0)
  mult <- vapply(tpval, function(v) sum(abs(tpval - v) < 1e-6), integer(1))
  tol <- 1e-5
  out <- vapply(uu, function(ui) {
    if (length(tp)) {
      j <- which.min(abs(ui - tp))
      if (abs(ui - tp[j]) < tol) {
        return(-2 * mult[j] * sign(polyval(cfsDD, tp[j])))
      }
    }
    y <- polyval(cfs, ui)
    r <- legendre_realroots(cfs, y)
    fY <- sum(1 / abs(polyval(cfsD, r)))
    fY * polyval(cfsD, ui)
  }, numeric(1))
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

# Breakpoints: 0, 1, and every non-smooth point of T. Similar in method to
# udplegendre's method, but the set of critical values whose other pre-images
# also need finding is bigger here: not just the turning-point values, but
# also g(0) and g(1) -- unlike a pure L_j, a general combination can have
# some OTHER u where g(u) = g(0) or g(u) = g(1) exactly (an ordinary point of
# g, but T inherits F's kink at that shared value all the same, exactly as
# it does at u = 0/1 themselves). If tp is empty, g is injective (no turning
# point means no repeated value by Rolle's theorem), so no cross points of
# either kind are possible and the early return is still exact.
setMethod("udpbreaks", "udplegendrebex", function(x) {
  cfs <- x@cfs
  cfsD <- x@cfsD
  tp <- legendre_turnpoints(cfsD)
  if (!length(tp)) {
    return(c(0, 1))
  }
  tol <- 5e-5
  yv <- sort(c(polyval(cfs, tp), polyval(cfs, c(0, 1))))
  yv <- yv[c(TRUE, diff(yv) > 1e-7)]
  base <- c(0, 1, tp)
  cross <- unlist(lapply(yv, function(y) {
    r <- legendre_realroots(cfs, y)
    r[vapply(r, function(ri) all(abs(ri - base) > tol), logical(1))]
  }))
  pts <- sort(c(base, cross))
  pts[c(TRUE, diff(pts) > tol)]
})

#' @describeIn pcoincide Integrate `sum_j p_j(v)^2` as in the default method,
#'   but split the range at the images of the turning points of `g`, where the
#'   integrand has a corner, so each piece is smooth.
#' @export
setMethod("pcoincide", "udplegendrebex", function(x) {
  tp <- legendre_turnpoints(x@cfsD)
  breaks <- pmin(pmax(x@Tfun(tp), 0), 1)
  integrate_collision(x, breaks)
})

#' Plot method for the udplegendrebex class
#'
#' Draws the graph of the udplegendrebex transformation over thin gridlines
#' marking its A-partition and T-partition: vertical lines at the points of
#' `udpbreaks()` (the A-partition, `0` and `1` excluded as redundant with the
#' plot's own border), the coarsest partition of `[0, 1]` into intervals on
#' which `T` is continuously differentiable; horizontal lines at the distinct
#' values `T` takes at those points (the T-partition), `0` and `1` again
#' excluded. A point with a one-sided vertical tangent but no fold reaches
#' the same `T` value as the turning point whose critical value it shares, so
#' it contributes a vertical line without adding a new horizontal one; `u = 0`
#' and `u = 1` are always A-partition members but need not be trivial
#' T-partition values, so they are handled separately in each role.
#'
#' @param x an object of class \linkS4class{udplegendrebex}.
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
#' plot(udplegendrebex(c(0.44, -0.35, 0.56, -0.21, 0.36)), embellish = "colour")
#' plot(udplegendrebex(c(0, 0.6, 0, -0.4)), embellish = "bw")
setMethod("plot", c(x = "udplegendrebex", y = "missing"),
  function(x, n = 500L, xlab = "u", ylab = "T(u)",
           embellish = c("none", "colour", "bw"), ...) {
    emb <- plot_embellish(embellish)
    b <- udpbreaks(x)
    u <- sort(unique(c(seq(0, 1, length.out = n), b)))
    # pty = "s" makes the plotting region square in physical inches; without
    # it, asp = 1 stretches one axis's displayed range past [0, 1] on any
    # device or panel that isn't already exactly square.
    op <- graphics::par(pty = "s")
    on.exit(graphics::par(op))
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
