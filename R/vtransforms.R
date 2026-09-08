#' Class of v-transforms
#'
#' This is the class of v-transforms. It contains the \linkS4class{vtransformi} subclass consisting of v-transforms
#' with an analytical expression for the inverse.
#'
#' @slot name a name for the v-transform of class character.
#' @slot vtrans function to evaluate the v-transform.
#' @slot pars vector containing the named parameters of the v-transform.
#' @slot gradient function to evaluate the gradient of the v-transform.
#'
#' @include udp-package.R
#' @export
#'
#' @examples
#' v2p(delta = 0.5, kappa = 1.2)
setClass("vtransform", contains = "udp", slots = list(
  name = "character", vtrans = "function", pars = "numeric",
  gradient = "function"
))

#' Class of invertible v-transforms
#'
#' This class inherits from the \linkS4class{vtransform} class and contains v-transforms
#' with an analytical expression for the inverse.
#'
#' @slot name a name for the v-transform of class character.
#' @slot vtrans function to evaluate the v-transform.
#' @slot pars vector containing the named parameters of the v-transform.
#' @slot gradient function to evaluate the gradient of the v-transform.
#' @slot inverse function to evaluate the inverse of the v-transform.
#'
#' @export
#'
#' @examples
#' vlinear(delta = 0.55)
setClass("vtransformi", contains = "vtransform", slots = list(
  name = "character", vtrans = "function",
  pars = "numeric", gradient = "function", inverse = "function"
))

# Validate a fulcrum: every v-transform must have delta strictly inside (0, 1).
check_delta <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
      delta <= 0 || delta >= 1) {
    stop("'delta' must be a single number in (0, 1).", call. = FALSE)
  }
}

#' Constructor function for symmetric v-transform
#'
#' @return An object of class \linkS4class{vtransformi}.
#' @export
#'
#' @examples
#' vsymmetric()
vsymmetric <- function() {
  new("vtransformi", name = "vsymmetric", vtrans = function(u) {
    abs(2 * u - 1)
  }, gradient = function(u) {
    slope <- rep(-2, length(u))
    slope[u > 0.5] <- 2
    slope
  }, inverse = function(v) {
    (1 - v) / 2
  })
}

#' Constructor function for linear v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#'
#' @return An object of class \linkS4class{vtransformi}.
#' @export
#'
#' @examples
#' vlinear(delta = 0.45)
vlinear <- function(delta = 0.5) {
  check_delta(delta)
  new("vtransformi", name = "vlinear", vtrans = function(u, delta) {
    abs(u / delta - 1) * ((delta / (1 - delta))^(u > delta))
  }, pars = c(delta = delta), gradient = function(u, delta) {
    slope <- rep(-1 / delta, length(u))
    slope[u > delta] <- 1 / (1 - delta)
    slope
  }, inverse = function(v, delta) {
    delta * (1 - v)
  })
}

#' Constructor function for 2-parameter v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{vtransform}.
#' @export
#'
#' @examples
#' v2p(delta = 0.45, kappa = 1.2)
v2p <- function(delta = 0.5, kappa = 1) {
  check_delta(delta)
  new("vtransform", name = "v2p", vtrans = function(u, delta, kappa) {
    ifelse(u <= delta, 1 - u - (1 - delta) * exp(-kappa * log(delta / u)), u - delta *
             exp(-(-log((1 - u) / (1 - delta)) / kappa)))
  }, pars = c(delta = delta, kappa = kappa),
  gradient = function(u, delta, kappa) {
    slope <- rep(-1, length(u))
    arg1 <- log(delta / u[u <= delta])
    arg2 <- log((1 - delta) / (1 - u[u > delta]))
    slope[u <= delta] <- -1 - (1 - delta) * exp(-kappa * arg1) * kappa / u[u <= delta]
    slope[u > delta] <- 1 + delta * exp(-arg2 / kappa) / (kappa * (1 - u[u > delta]))
    if (length(u[u == 0]) > 0) {
      if (kappa < 1) {
        val <- -Inf
      }
      if (kappa > 1) {
        val <- -1
      }
      if (kappa == 1) {
        val <- -1 / delta
      }
      slope[(u == 0)] <- val
    }
    if (length(u[u == 1]) > 0) {
      if (kappa > 1) {
        val2 <- Inf
      }
      if (kappa < 1) {
        val2 <- 1
      }
      if (kappa == 1) {
        val2 <- 1 / (1 - delta)
      }
      slope[(u == 1)] <- val2
    }
    slope
  })
}

#' Constructor function for 2-parameter beta v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{vtransform}.
#' @export
#'
#' @examples
#' v2b(delta = 0.45, kappa = 1.2)
v2b <- function(delta = 0.5, kappa = 1) {
  check_delta(delta)
  new("vtransform", name = "v2b", vtrans = function(u, delta, kappa) {
    suppressWarnings(ifelse(u <= delta, 1 - u - (1 - delta) * pbeta(
      u / delta, kappa,
      1 / kappa
    ), u - delta * qbeta((1 - u) / (1 - delta), kappa, 1 / kappa)))
  }, pars = c(delta = delta, kappa = kappa),
  gradient = function(u, delta, kappa) {
    slope <- rep(-1, length(u))
    slope[u <= delta] <- -1 - (1 - delta) * dbeta(u[u <= delta] / delta, kappa, 1 / kappa) / delta
    slope[u > delta] <- 1 + delta / (dbeta(qbeta((1 - u[u > delta]) / (1 - delta), kappa, 1 / kappa), kappa, 1 / kappa) * (1 - delta))
    slope
  })
}

#' Constructor function for 3-parameter v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#' @param xi additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{vtransform}.
#' @export
#'
#' @examples
#' v3p(delta = 0.45, kappa = 0.8, xi = 1.1)
v3p <- function(delta = 0.5, kappa = 1, xi = 1) {
  check_delta(delta)
  new("vtransform", name = "v3p", vtrans = function(u, delta, kappa, xi) {
    ifelse(u <= delta, 1 - u - (1 - delta) * exp(-kappa * (log(delta / u))^xi), u -
             delta * exp(-(-log((1 - u) / (1 - delta)) / kappa)^(1 / xi)))
  }, pars = c(delta = delta, kappa = kappa, xi = xi),
  gradient = function(u, delta, kappa, xi) {
    slope <- rep(-1, length(u))
    arg1 <- log(delta / u[u <= delta])
    arg2 <- log((1 - delta) / (1 - u[u > delta]))
    slope[u <= delta] <- -1 - (1 - delta) * exp(-kappa * arg1^xi) * arg1^(xi -
                                                                            1) * (xi * kappa / u[u <= delta])
    slope[u > delta] <- 1 + delta * exp(-(arg2 / kappa)^(1 / xi)) *
      (arg2 / kappa)^(1 / xi - 1) / (xi * kappa * (1 - u[u > delta]))
    if (length(u[u == 0]) > 0) {
      if ((xi < 1) || ((xi == 1) && (kappa < 1))) {
        val <- -Inf
      }
      if ((xi > 1) || ((xi == 1) && (kappa > 1))) {
        val <- -1
      }
      if ((xi == 1) && (kappa == 1)) {
        val <- -1 / delta
      }
      slope[(u == 0)] <- val
    }
    if (length(u[u == 1]) > 0) {
      if ((xi > 1) || ((xi == 1) && (kappa < 1))) {
        val2 <- Inf
      }
      if ((xi < 1) || ((xi == 1) && (kappa > 1))) {
        val2 <- 1
      }
      if ((xi == 1) && (kappa == 1)) {
        val2 <- 1 / (1 - delta)
      }
      slope[(u == 1)] <- val2
    }
    slope
  })
}

#' Constructor function for 3-parameter beta v-transform
#'
#' @param delta a value in (0, 1) specifying the fulcrum of the v-transform.
#' @param kappa additional positive parameter of v-transform.
#' @param xi additional positive parameter of v-transform.
#'
#' @return An object of class \linkS4class{vtransform}.
#' @export
#'
#' @examples
#' v3b(delta = 0.45, kappa = 1.2, xi = 1.2)
v3b <- function(delta = 0.5, kappa = 1, xi = 1) {
  check_delta(delta)
  new("vtransform", name = "v3b", vtrans = function(u, delta, kappa, xi) {
    suppressWarnings(ifelse(u <= delta, 1 - u - (1 - delta) * pbeta(
      u / delta, kappa,
      xi
    ), u - delta * qbeta((1 - u) / (1 - delta), kappa, xi)))
  }, pars = c(delta = delta, kappa = kappa, xi = xi),
  gradient = function(u, delta, kappa, xi) {
    slope <- rep(-1, length(u))
    slope[u <= delta] <- -1 - (1 - delta) * dbeta(u[u <= delta] / delta, kappa, xi) / delta
    slope[u > delta] <- 1 + delta / (dbeta(qbeta((1 - u[u > delta]) / (1 - delta), kappa, xi), kappa, xi) * (1 - delta))
    slope
  })
}

#' @describeIn udptrans Evaluate a v-transform.
#' @export
setMethod("udptrans", "vtransform", function(x, u) {
  do.call(x@vtrans, append(x@pars, list(u = u)))
})

#' Calculate gradient of v-transform
#'
#' @param x an object of class \linkS4class{vtransform}.
#' @param u a vector, matrix or time series with values in `[0, 1]`.
#'
#' @return An object shaped like `u` giving the gradient of the v-transform.
#' @export
#'
#' @examples
#' vgradient(vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
vgradient <- function(x, u) {
  g <- do.call(x@gradient, append(x@pars, list(u = u)))
  if (!is.null(attributes(u))) {
    attributes(g) <- attributes(u)
  }
  g
}

#' Calculate the lower-branch inverse of a v-transform
#'
#' Returns the pre-image at or below the fulcrum: the value `u` in
#' `[0, delta]` with `udptrans(x, u)` equal to `v`.
#'
#' For a \linkS4class{vtransformi} object the analytic inverse stored in the
#' `inverse` slot is used and `method`, `tol` and `ngrid` are ignored.
#' Otherwise the inverse is computed numerically, either with
#'
#' * `method = "newton"` (the default): a vectorised Newton iteration,
#'   safeguarded by bisection, that uses the analytic gradient of the
#'   v-transform. Accurate to roughly `tol` and typically several times
#'   faster than element-wise root finding.
#' * `method = "spline"`: monotone cubic interpolation of the v-transform
#'   evaluated on an equally spaced grid of `ngrid` points on `[0, delta]`.
#'   Much faster again for long `v`, but the accuracy is limited by the grid
#'   spacing and degrades for extreme parameter values (for example a small
#'   `kappa`, where the v-transform has infinite slope at `0`).
#'
#' @param x an object of class \linkS4class{vtransform}.
#' @param v a vector, matrix or time series with values in `[0, 1]`.
#' @param method inversion method for non-invertible v-transforms, either
#' `"newton"` or `"spline"`. Ignored for invertible v-transforms.
#' @param tol convergence tolerance for `method = "newton"`.
#' @param ngrid number of grid points (at least 2) for `method = "spline"`.
#'
#' @return An object shaped like `v` with values in `[0, delta]`. Positions
#' where `v` is `NA` or otherwise non-finite are returned as `NA`.
#' @export
#'
#' @examples
#' vinverse(vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
#' vinverse(v2p(delta = 0.4, kappa = 1.3), seq(0.1, 0.9, by = 0.2))
#' vinverse(v2p(delta = 0.4, kappa = 1.3), seq(0.1, 0.9, by = 0.2), method = "spline")
vinverse <- function(x, v, method = c("newton", "spline"),
                     tol = .Machine$double.eps^0.5, ngrid = 1000L) {
  method <- match.arg(method)

  if (is(x, "vtransformi")) {
    return(do.call(x@inverse, append(x@pars, list(v = v))))
  }

  delta <- unname(x@pars["delta"])
  if (is.na(delta)) {
    stop("'x' has no 'delta' parameter; cannot invert numerically.")
  }
  parlist <- as.list(x@pars)
  vf <- function(u) do.call(x@vtrans, c(list(u = u), parlist))

  vv <- as.numeric(v)
  out <- rep(NA_real_, length(vv))
  finite <- is.finite(vv)
  out[finite & vv <= 0] <- delta
  out[finite & vv >= 1] <- 0
  todo <- finite & vv > 0 & vv < 1

  if (any(todo)) {
    vt <- vv[todo]
    if (method == "spline") {
      if (length(ngrid) != 1L || !is.finite(ngrid) || ngrid < 2) {
        stop("'ngrid' must be a single number of at least 2.")
      }
      ug <- seq(0, delta, length.out = ngrid)
      out[todo] <- stats::splinefun(rev(vf(ug)), rev(ug), method = "monoH.FC")(vt)
    } else {
      vg <- function(u) do.call(x@gradient, c(list(u = u), parlist))
      lo <- rep(0, length(vt))
      hi <- rep(delta, length(vt))
      u <- rep(delta / 2, length(vt))
      f <- vf(u) - vt
      g <- vg(u)
      dx <- dxold <- rep(delta, length(vt))
      for (i in seq_len(100L)) {
        # phi(u) = udptrans(x, u) - v is strictly decreasing on [0, delta]
        left <- f > 0
        lo[left] <- u[left]
        hi[!left] <- u[!left]

        step_ok <- is.finite(g) & g != 0
        cand <- u - f / g
        slow <- abs(2 * f) > abs(dxold * g)
        bisect <- !step_ok | slow | cand <= lo | cand >= hi
        cand[bisect] <- 0.5 * (lo[bisect] + hi[bisect])

        dxold <- dx
        dx <- cand - u
        u <- cand
        if (max(abs(dx)) < tol) break

        f <- vf(u) - vt
        g <- vg(u)
      }
      if (max(abs(dx)) >= tol) {
        warning("vinverse(): Newton iteration did not reach 'tol' in 100 steps.")
      }
      out[todo] <- u
    }
  }

  if (!is.null(attributes(v))) {
    attributes(out) <- attributes(v)
  }
  out
}

#' Calculate conditional down probability of v-transform
#'
#' @param x an object of class \linkS4class{vtransform}.
#' @param v a vector or time series with values in `[0, 1]`.
#' @param tol convergence tolerance passed to [vinverse()].
#' @param ... further arguments passed to [vinverse()], such as `method`.
#'
#' @return A vector or time series of values of gradient.
#' @export
#'
#' @examples
#' vdownprob(v2p(delta = 0.55, kappa = 1.2), c(0, 0.25, 0.5, 0.75, 1))
vdownprob <- function(x, v, tol = .Machine$double.eps^0.5, ...) {
  -1 / vgradient(x, vinverse(x, v, tol = tol, ...))
}

#' @describeIn udpsi Stochastic inverse of a v-transform. Accepts `tol` and
#'   further arguments of [vinverse()] (such as `method`).
#' @export
setMethod("udpsi", "vtransform", function(x, v, Z = runif(length(v)),
                                          tol = .Machine$double.eps^0.5, ...) {
  if (length(Z) != length(v)) {
    stop("'Z' must have the same length as 'v'.")
  }
  vinv <- vinverse(x, v, tol = tol, ...)
  pdown <- -1 / vgradient(x, vinv)
  # the two pre-images of v satisfy u2 = u1 + v, so v + vinv is the upper one
  output <- ifelse(Z <= pdown, vinv, v + vinv)
  if (!(is.null(attributes(v)))) {
    attributes(output) <- attributes(v)
  }
  output
})

#' Plot method for vtransform class
#'
#' Plots the v-transform as well as its gradient or inverse. Can also plot the
#' conditional probability that a series PIT falls below the fulcrum for a
#' given volatility PIT value v.
#'
#' @param x an object of class \linkS4class{vtransform}.
#' @param type type of plot: 'transform' for plot of transform, 'inverse' for plot of inverse,
#' 'gradient' for plot of gradient or 'pdown' for plot of conditional probability.
#' @param shading logical variable specifying whether inadmissible zone for v-transform
#' should be shaded
#' @param npoints number of plotting points along x-axis.
#' @param lower the lower x-axis value for plotting.
#' @param upper the upper x-axis value for plotting
#'
#' @return No return value, generates plot.
#' @export
#'
#'
#' @examples
#' plot(vsymmetric())
#' plot(v2p(delta = 0.45, kappa = 0.8), type = "inverse")
#' plot(v2p(delta = 0.45, kappa = 0.8), type = "gradient")
setMethod("plot", c(x = "vtransform", y = "missing"), function(x, type = "transform",
                                                               shading = TRUE, npoints = 200, lower = 0, upper = 1) {
  delta <- ifelse(is.element("delta", names(x@pars)), x@pars["delta"], 0.5)
  switch(type, inverse = {
    vvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    plot(vvals, vinverse(x, vvals), xlab = "v", ylab = "Vinv(v)", type = "l")
  }, gradient = {
    uvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    plot(uvals, vgradient(x, uvals), xlab = "u", ylab = "Vprime(u)", type = "l")
  }, pdown = {
    vvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    plot(vvals, vdownprob(x, vvals), xlab = "v", ylab = "Delta(v)", type = "l")
  }, transform = {
    uvals <- seq(from = max(lower, 0), to = min(upper, 1), length = npoints)
    if ((delta > lower) & (delta < upper)) uvals <- sort(c(uvals, delta))
    plot(uvals, udptrans(x, uvals), xlab = "u", ylab = "V(u)", type = "l")
    if (shading) {
      # colchoice = 'gray97'
      colchoice <- "gray90"
      polygon(c(0, 0, delta), c(delta, 0, 0), col = colchoice, border = NA)
      polygon(c(delta, 1, 1), c(0, 0, 1 - delta), col = colchoice, border = NA)
      polygon(c(0, delta, delta), c(1, 1, 1 - delta), col = colchoice, border = NA)
      polygon(c(delta, delta, 1), c(delta, 1, 1), col = colchoice, border = NA)
    }
  }, stop("Not a plot method for v-transform."))
})

#' Compute coincidence probability for v-transform
#'
#' Computes the probability that if we v-transform a uniform
#' random variable and then stochastically invert the
#' v-transform, we get back to the original value.
#'
#' @param x an object of class \linkS4class{vtransform}.
#'
#' @return The probability of coincidence.
#' @export
#'
#' @examples
#' pcoincide(vlinear(delta = 0.4))
#' pcoincide(v3p(delta = 0.45, kappa = 0.5, xi = 1.3))
pcoincide <- function(x) {
  if (x@name == "vsymmetric") {
    return(0.5)
  }
  delta <- unname(x@pars["delta"])
  # the linear v-transform has constant down-probability delta, so varDelta = 0
  if (x@name == "vlinear") {
    return(delta^2 + (1 - delta)^2)
  }
  integrand <- function(v) (vdownprob(x, v) - delta)^2
  varDelta <- integrate(integrand, 0, 1)$value
  unname(delta^2 + (1 - delta)^2 + 2 * varDelta)
}

#' @describeIn vtransform Show method for vtransform class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", "vtransform", function(object) {
  cat("name: ", object@name, "\n", "parameters: ",
      "\n",
      sep = ""
  )
  print(object@pars)
})

#' @describeIn vtransform Coef method for vtransform class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", "vtransform", function(object) {
  object@pars
})
