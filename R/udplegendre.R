# udp transformations built from shifted Legendre polynomials.
#
# For a shifted Legendre polynomial L_j on [0, 1] and U uniform on [0, 1], let
# F_j be the distribution function of L_j(U). The transformation is the
# composition T(u) = F_j(L_j(u)), which is uniform-distribution-preserving:
# F_j is strictly increasing on the range of L_j, so T has exactly the (up to
# j) pre-images of L_j, and the stochastic-inverse weights 1 / |T'| reduce to
# 1 / |L_j'| after normalisation (the common F_j' factor cancels).

# Ascending monomial coefficients of the shifted Legendre polynomial of the
# given degree on [0, 1], from the three-term recurrence
#   L_0 = 1,  L_1 = 2x - 1,
#   (n + 1) L_{n+1} = (2n + 1)(2x - 1) L_n - n L_{n-1}.
# Coefficients are exact integers well past the usable degree range; for degree
# beyond roughly 12 the alternating monomial coefficients make polyroot() and
# Horner evaluation lose accuracy.
slegendre_coef <- function(degree) {
  prev <- 1 # L_0
  if (degree == 0L) {
    return(prev)
  }
  cur <- c(-1, 2) # L_1
  for (n in seq_len(degree - 1L)) {
    # (2x - 1) * cur in ascending coefficients: 2x * cur raises every power,
    # -1 * cur leaves it
    shifted <- c(0, 2 * cur) + c(-cur, 0)
    prev <- c(prev, numeric(length(shifted) - length(prev)))
    nxt <- ((2 * n + 1) * shifted - n * prev) / (n + 1)
    prev <- cur
    cur <- nxt
  }
  cur
}

# Ascending coefficients of the derivative of a polynomial given by ascending
# coefficients: d/dx sum a_k x^k = sum k a_k x^{k-1}.
poly_deriv_coef <- function(coef) {
  if (length(coef) <= 1L) {
    return(0)
  }
  coef[-1L] * seq_len(length(coef) - 1L)
}

# Horner evaluation of a polynomial given by ascending coefficients.
polyval <- function(coef, x) {
  y <- numeric(length(x))
  for (a in rev(coef)) {
    y <- y * x + a
  }
  y
}

# Real roots in [0, 1] of L_j(u) = y, where `coef` are the ascending
# coefficients of L_j. Candidates from polyroot() are kept when they sit near
# the real axis and near [0, 1], then verified by residual against the local
# slope; near-coincident roots (y at a turning point) are merged.
legendre_realroots <- function(coef, y, tol = 1e-7) {
  cc <- coef
  cc[1] <- cc[1] - y
  z <- polyroot(cc)
  re <- Re(z)
  cand <- re[abs(Im(z)) < 1e-6 * pmax(1, abs(re)) & re > -tol & re < 1 + tol]
  if (!length(cand)) {
    return(numeric(0))
  }
  cand <- pmin(pmax(cand, 0), 1)
  resid <- abs(polyval(coef, cand) - y)
  slope <- pmax(abs(polyval(poly_deriv_coef(coef), cand)), 1)
  cand <- sort(cand[resid < 1e-6 * slope + 1e-9])
  if (length(cand) > 1L) {
    cand <- cand[c(TRUE, diff(cand) > tol)]
  }
  cand
}

# Interior turning points of L_j (roots of L_j' in (0, 1)); all j - 1 of them
# are real.
legendre_turnpoints <- function(coefD) {
  if (length(coefD) < 2L) {
    return(numeric(0))
  }
  z <- polyroot(coefD)
  tp <- Re(z)[abs(Im(z)) < 1e-6]
  sort(tp[tp > 0 & tp < 1])
}

# Lower end of the range of L_j on [0, 1]: -1 for odd degree (attained at
# u = 0), the smallest interior turning-point value for even degree.
legendre_lbound <- function(coef, coefD, degree) {
  if (degree %% 2L == 1L) {
    return(-1)
  }
  min(polyval(coef, c(legendre_turnpoints(coefD), 0, 1)))
}

# Exact F_j(y) = |{u in [0, 1] : L_j(u) <= y}| for a single y. Between
# consecutive roots of L_j(u) - y the polynomial keeps one side of y, so the
# measure is the total length of the sub-intervals whose midpoint sits at or
# below y. Robust to a dropped near-turning-point root (it bounds a vanishing
# interval) and needs no assumption on the number of roots.
legendre_measure <- function(coef, y, lbound) {
  if (y <= lbound) {
    return(0)
  }
  if (y >= 1) {
    return(1)
  }
  breaks <- c(0, legendre_realroots(coef, y), 1)
  mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  sum(diff(breaks)[polyval(coef, mids) <= y])
}

#' Class of shifted-Legendre udp transformations
#'
#' A shifted-Legendre udp transformation is `T(u) = F_j(L_j(u))`, the
#' composition of the shifted Legendre polynomial `L_j` of the stated degree
#' `j` on `[0, 1]` with the distribution function `F_j` of `L_j(U)` for `U`
#' uniform. It is a degree-`j`, uniform-distribution-preserving map of
#' `[0, 1]`: for `U` uniform, `udptrans(x, U)` is again uniform.
#'
#' `F_j` and its inverse have no closed form for `degree` above 2, so the
#' constructor evaluates `F_j` exactly on a grid once (from the real roots of
#' `L_j(u) = y`) and stores spline interpolations, `Tfun` for `T` itself and
#' `Qfun` for the quantile `F_j^{-1}`. `F_j` has a square-root branch point at
#' every interior extreme value of `L_j`; the splines are built panel by panel
#' between those values under the substitution that regularises them, so the
#' interpolation error stays near `1e-4` across the usable degree range.
#' Pre-image finding in [udpinverse()] stays exact (polynomial root finding);
#' only the `F_j` reparametrisation is interpolated.
#'
#' @slot degree integer; the degree of the transformation (at least 1).
#' @slot cfs,cfsD ascending monomial coefficients of `L_j` and of its
#'   derivative.
#' @slot lbound lower end of the range of `L_j` on `[0, 1]`.
#' @slot Tfun,Qfun functions evaluating `T` and `F_j^{-1}`.
#'
#' @seealso [udplegendre()] to construct one, [udptrans()] to evaluate it.
#' @include udp-package.R
#' @export
setClass("udplegendre",
  contains = "udp",
  slots = list(
    degree = "integer", cfs = "numeric", cfsD = "numeric",
    lbound = "numeric", Tfun = "function", Qfun = "function"
  )
)

#' Construct a shifted-Legendre udp transformation
#'
#' @param degree a single positive integer, the degree of the transformation.
#' @param ngrid number of grid points (at least 3) at which `F_j` is evaluated
#'   exactly before interpolation; ignored for `degree` 1 and 2, which are
#'   closed form.
#'
#' @return An object of class \linkS4class{udplegendre}.
#' @export
#'
#' @examples
#' udplegendre(3)
#' plot(function(u) udptrans(udplegendre(4), u), xlab = "u", ylab = "T(u)")
udplegendre <- function(degree, ngrid = 257L) {
  if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) ||
    degree < 1 || degree != round(degree)) {
    stop("'degree' must be a single positive integer.", call. = FALSE)
  }
  degree <- as.integer(degree)
  if (degree > 12L) {
    warning(
      "udplegendre(): degree > 12 is numerically unreliable (polyroot on ",
      "large alternating coefficients).",
      call. = FALSE
    )
  }
  cfs <- slegendre_coef(degree)
  cfsD <- poly_deriv_coef(cfs)

  if (degree == 1L) {
    lbound <- -1
    Tfun <- function(u) u
    Qfun <- function(v) 2 * v - 1
  } else if (degree == 2L) {
    lbound <- -0.5
    Tfun <- function(u) abs(2 * u - 1)
    Qfun <- function(v) (3 * v^2 - 1) / 2
  } else {
    if (!is.numeric(ngrid) || length(ngrid) != 1L || is.na(ngrid) || ngrid < 3) {
      stop("'ngrid' must be a single number of at least 3.", call. = FALSE)
    }
    ngrid <- as.integer(ngrid)
    lbound <- legendre_lbound(cfs, cfsD, degree)
    # interior extreme values of L_j, symmetric pairs merged
    yv <- sort(polyval(cfs, legendre_turnpoints(cfsD)))
    yv <- yv[c(TRUE, diff(yv) > 1e-7)]
    ye <- yv[yv > lbound + 1e-7 & yv < 1 - 1e-7]
    ypanels <- c(lbound, ye, 1)
    npan <- length(ypanels) - 1L
    per <- max(40L, ceiling(ngrid / npan))

    # F_j has a square-root branch point at every interior extreme value of
    # L_j. On panel p = [a, b] between consecutive such values the substitution
    # y = mid - half * cos(theta) regularises the sqrt(|y - a|) and
    # sqrt(|b - y|) behaviour at the two ends, so F_j is smooth in theta and a
    # modest monotone spline (and its inverse) capture it. udptrans() then
    # evaluates F_j(L_j(u)); udpinverse() evaluates F_j^{-1}.
    panel <- lapply(seq_len(npan), function(p) {
      a <- ypanels[p]
      b <- ypanels[p + 1L]
      mid <- (a + b) / 2
      half <- (b - a) / 2
      th <- seq(0, pi, length.out = per)
      yy <- mid - half * cos(th)
      yy[c(1L, per)] <- c(a, b)
      fv <- cummax(vapply(yy, legendre_measure, numeric(1),
        coef = cfs, lbound = lbound
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
      pmin(pmax(panel_apply(pmin(pmax(y, lbound), 1), ypanels, function(pl, yy) {
        pl$Fofth(acos(pmin(pmax((pl$mid - yy) / pl$half, -1), 1)))
      }), 0), 1)
    }
    Qfun <- function(v) {
      pmin(pmax(panel_apply(pmin(pmax(v, 0), 1), vpanels, function(pl, vv) {
        pl$mid - pl$half * cos(pmin(pmax(pl$thofF(vv), 0), pi))
      }), lbound), 1)
    }
    Tfun <- function(u) Ffun(polyval(cfs, pmin(pmax(u, 0), 1)))
  }

  new("udplegendre",
    degree = degree, cfs = cfs, cfsD = cfsD, lbound = lbound,
    Tfun = Tfun, Qfun = Qfun
  )
}

#' @describeIn udptrans Evaluate a shifted-Legendre udp transformation.
#' @export
setMethod("udptrans", "udplegendre", function(x, u) {
  out <- pmin(pmax(x@Tfun(pmin(pmax(as.numeric(u), 0), 1)), 0), 1)
  if (!is.null(attributes(u))) {
    attributes(out) <- attributes(u)
  }
  out
})

#' @describeIn udpinverse Pre-images of a shifted-Legendre udp transformation:
#'   a matrix with `degree` columns holding, for each `v`, the roots in
#'   `[0, 1]` of `L_j(u) = F_j^{-1}(v)`, sorted ascending and left-packed with
#'   trailing `NA`. With `prob = TRUE` the `"prob"` attribute weights each root
#'   by `1 / |L_j'(u)|`, normalised over the row.
#' @export
setMethod("udpinverse", "udplegendre", function(x, v, prob = FALSE, ...) {
  vv <- as.numeric(v)
  if (anyNA(vv) || any(vv < 0 | vv > 1)) {
    stop("every element of 'v' must be in [0, 1].", call. = FALSE)
  }
  k <- x@degree
  y <- pmin(pmax(x@Qfun(vv), x@lbound), 1)
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

#' @describeIn pcoincide Integrate `sum_j p_j(v)^2` as in the default method,
#'   but split the range at the images of the turning points of `L_j`, where
#'   the integrand has a corner, so each piece is smooth. Accuracy is bounded
#'   by the `F_j^{-1}` spline, like the rest of the class.
#' @export
setMethod("pcoincide", "udplegendre", function(x) {
  tp <- legendre_turnpoints(x@cfsD)
  breaks <- pmin(pmax(x@Tfun(tp), 0), 1)
  integrate_collision(x, breaks)
})

#' Plot method for the udplegendre class
#'
#' Draws the graph of the shifted-Legendre udp transformation over thin red
#' gridlines: vertical at the turning points of `L_j` (which are also the
#' turning points of `T`) and horizontal at the transformed turning-point
#' values `T(tp)`.
#'
#' @param x an object of class \linkS4class{udplegendre}.
#' @param n number of points at which to evaluate the transformation.
#' @param xlab,ylab axis labels.
#' @param ... further graphical parameters passed to [graphics::plot()].
#'
#' @return No return value, generates a plot.
#' @export
#'
#' @examples
#' plot(udplegendre(4))
#' plot(udplegendre(7))
setMethod("plot", c(x = "udplegendre", y = "missing"),
  function(x, n = 500L, xlab = "u", ylab = "T(u)", ...) {
    tp <- legendre_turnpoints(x@cfsD)
    u <- sort(unique(c(seq(0, 1, length.out = n), tp)))
    plot(NA,
      xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
      xlab = xlab, ylab = ylab, ...
    )
    abline(v = tp, h = udptrans(x, tp), col = "red", lwd = 0.5)
    lines(u, udptrans(x, u), lwd = 2)
  }
)
