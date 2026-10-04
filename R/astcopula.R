# The absolute spherical t copula (ast): the copula of (|X1|, |X2|) for a
# bivariate t distribution with nu degrees of freedom and correlation 0.
# Equivalently, the t copula with rho = 0 transformed on both margins by the
# symmetric v-transform T(u) = |2u - 1|; applying the stochastic inverse of
# T to both margins (independent randomizers) gives that t copula back.
#
# With p = (1 + u) / 2, q = (1 + v) / 2, a = t_nu^-1(p) and b = t_nu^-1(q),
# the symmetry of the t copula with rho = 0 under a sign change of either
# coordinate gives (Dias, Han and McNeil, Supplementary Material, "The
# absolute spherical t copula"):
#
#   density  c(u, v)  = c_t(p, q)
#   CDF      C(u, v)  = 4 C_t(p, q) - u - v - 1
#   h        dC/du    = 2 t_{nu+1}( b sqrt((nu + 1) / (nu + a^2)) ) - 1
#
# (the printed h-function there omits the square on a). Given X1 = a, X2 is
# t_{nu+1} scaled by sqrt((nu + a^2) / (nu + 1)), so h and its inverse need
# only the univariate t. The CDF is the integral of h, since rvinecopulib's t
# copula CDF needs nu >= 2 and copula's integer nu.
#
# All quantiles are taken in the upper tail, t_nu^-1((1 + u) / 2) =
# qt((1 - u) / 2, nu, lower.tail = FALSE), so u close to 1 keeps its accuracy.

# Upper-tail quantile a(u) = t_nu^-1((1 + u) / 2), Inf at u = 1.
ast_a <- function(u, nu) stats::qt((1 - u) / 2, nu, lower.tail = FALSE)

# 2 t_{nu+1}(b sqrt((nu + 1) / (nu + a^2))) - 1, with a, b >= 0 possibly Inf.
ast_h_ab <- function(a, b, nu) {
  z <- b * sqrt((nu + 1) / (nu + a^2))
  z[is.infinite(b)] <- Inf
  z[is.infinite(a) & !is.infinite(b)] <- 0
  1 - 2 * stats::pt(z, nu + 1, lower.tail = FALSE)
}

ast_density <- function(u, v, nu) {
  x <- ast_a(u, nu)
  y <- ast_a(v, nu)
  logc <- lgamma((nu + 2) / 2) - lgamma(nu / 2) - log(nu * pi) -
    (nu + 2) / 2 * log1p((x^2 + y^2) / nu) -
    stats::dt(x, nu, log = TRUE) - stats::dt(y, nu, log = TRUE)
  exp(logc)
}

# h(v | u) = P(V <= v | U = u), the derivative of the CDF in u; given_var = 2
# swaps the roles. Exact at v = 0 and v = 1.
ast_h <- function(given, x, nu) ast_h_ab(ast_a(given, nu), ast_a(x, nu), nu)

# The v with ast_h(given, v, nu) = w.
ast_hinv <- function(given, w, nu) {
  a <- ast_a(given, nu)
  b <- sqrt((nu + a^2) / (nu + 1)) * stats::qt((1 - w) / 2, nu + 1, lower.tail = FALSE)
  1 - 2 * stats::pt(b, nu, lower.tail = FALSE)
}

# C(u, v) by integrating h over the smaller argument (C is symmetric), with
# tanh-sinh quadrature: the integrand is smooth but its derivative blows up
# at s = 1.
ast_cdf <- function(u, v, nu, N = 41L) {
  lo <- pmin(u, v)
  hi <- pmax(u, v)
  b <- ast_a(hi, nu)
  out <- cell_quad(0, lo, N, function(s, idx) ast_h_ab(ast_a(s, nu), b[idx], nu))
  out[lo >= 1] <- hi[lo >= 1]
  pmin(pmax(out, 0), pmin(lo, hi))
}

ast_sample <- function(n, nu) {
  x <- abs(matrix(stats::rnorm(2L * n), n, 2L)) / sqrt(stats::rchisq(n, nu) / nu)
  1 - 2 * stats::pt(x, nu, lower.tail = FALSE)
}

#' The absolute spherical t copula
#'
#' The copula of `(|X1|, |X2|)`, where `(X1, X2)` has a bivariate t
#' distribution with `nu` degrees of freedom and correlation zero. It is the
#' t copula with correlation zero transformed on both margins by the
#' symmetric v-transform ([vsymmetric()]), and applying the stochastic inverse
#' of that transformation under independent randomizers gives the t copula
#' back (Dias, Han and McNeil, 2027). The family interpolates between
#' independence (`nu` to infinity) and comonotonicity (`nu` to 0), so
#' Kendall's tau ([astcopula_tau()]) ranges over `(0, 1)`. It has upper tail
#' dependence and is asymptotically independent in the lower tail.
#'
#' `astcopula()` creates the object. It can be used as the `basecopula` of a
#' \linkS4class{bsicopula}, but not in a \linkS4class{randsdvine}.
#' `dastcopula()`, `pastcopula()`, `hastcopula()` and `rastcopula()` give the
#' density, distribution function, h-functions and their inverses, and random
#' numbers. With `p = (1 + u) / 2` and `q = (1 + v) / 2`, the density is the
#' t copula density at `(p, q)` and the CDF is `4 C_t(p, q) - u - v - 1`.
#' The h-function is `2 t_{nu+1}(b sqrt((nu + 1) / (nu + a^2))) - 1`, where
#' `a` and `b` are the quantiles `t_nu^-1(p)` and `t_nu^-1(q)`, and its
#' inverse is also in closed form, so both need only the univariate t
#' distribution. The CDF is computed by integrating the h-function
#' (tanh-sinh quadrature), which is valid for every `nu > 0`; the t copula
#' CDF available in other packages restricts `nu`.
#'
#' @param nu degrees of freedom, a positive number.
#' @param u1,u2 numeric vectors with values in `[0, 1]`, recycled to a common
#'   length; `u1` may instead be a two-column matrix with `u2` omitted.
#' @param object an object created by `astcopula()`.
#' @param cond_var the conditioning variable of `hastcopula()`, `1` or `2`
#'   (the convention of [rvinecopulib::hbicop()]).
#' @param inverse logical; return the inverse h-function (the non-conditioned
#'   argument is then a probability level)?
#' @param n number of random pairs.
#'
#' @return `astcopula()` an object of class `astcopula`; the others numeric
#'   vectors, except `rastcopula()`, which gives an `n` by 2 matrix.
#' @references
#' Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#' D-vines for stochastic volatility. *Journal of Multivariate Analysis*,
#' **217**, 105695. \doi{10.1016/j.jmva.2026.105695}
#' @seealso [astcopula_tau()] and [astcopula_nu()] for Kendall's tau and the
#'   calibration of `nu` to it.
#' @include basecopula.R
#' @rdname astcopula
#' @aliases astcopula-class
#' @export
#' @examples
#' ac <- astcopula(3)
#' dastcopula(c(0.2, 0.9), c(0.3, 0.95), ac)
#' pastcopula(0.5, 0.5, ac)
#' h <- hastcopula(0.4, c(0.2, 0.6), ac)
#' hastcopula(0.4, h, ac, inverse = TRUE)
#' head(rastcopula(5, ac))
setClass("astcopula", slots = list(nu = "numeric"), validity = function(object) {
  if (length(object@nu) != 1L || is.na(object@nu) || !is.finite(object@nu) || object@nu <= 0) {
    "'nu' must be a single positive number."
  } else {
    TRUE
  }
})

#' @rdname astcopula
#' @export
astcopula <- function(nu) new("astcopula", nu = as.numeric(nu))

ast_args <- function(u1, u2, object) {
  if (!methods::is(object, "astcopula")) {
    stop("'object' must be an object of class 'astcopula'.", call. = FALSE)
  }
  if (is.null(u2)) {
    if (!(is.matrix(u1) || is.data.frame(u1)) || ncol(u1) != 2L) {
      stop("when 'u2' is omitted, 'u1' must be a two-column matrix.", call. = FALSE)
    }
    u1 <- as.matrix(u1)
    u2 <- u1[, 2L]
    u1 <- u1[, 1L]
  }
  u1 <- as.numeric(u1)
  u2 <- as.numeric(u2)
  n <- max(length(u1), length(u2))
  if (length(u1) == 1L) u1 <- rep(u1, n)
  if (length(u2) == 1L) u2 <- rep(u2, n)
  if (length(u1) != length(u2)) {
    stop("'u1' and 'u2' must have the same length (or one of them length 1).", call. = FALSE)
  }
  if (anyNA(u1) || anyNA(u2) || any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1)) {
    stop("'u1' and 'u2' must be numbers in [0, 1], without missing values.", call. = FALSE)
  }
  list(u1 = u1, u2 = u2, nu = object@nu)
}

#' @rdname astcopula
#' @export
dastcopula <- function(u1, u2 = NULL, object) {
  a <- ast_args(u1, u2, object)
  ast_density(a$u1, a$u2, a$nu)
}

#' @rdname astcopula
#' @export
pastcopula <- function(u1, u2 = NULL, object) {
  a <- ast_args(u1, u2, object)
  ast_cdf(a$u1, a$u2, a$nu)
}

#' @rdname astcopula
#' @export
hastcopula <- function(u1, u2 = NULL, object, cond_var = 1L, inverse = FALSE) {
  a <- ast_args(u1, u2, object)
  if (!(length(cond_var) == 1L && cond_var %in% 1:2)) {
    stop("'cond_var' must be 1 or 2.", call. = FALSE)
  }
  g <- if (cond_var == 1L) a$u1 else a$u2
  x <- if (cond_var == 1L) a$u2 else a$u1
  if (inverse) ast_hinv(g, x, a$nu) else ast_h(g, x, a$nu)
}

#' @rdname astcopula
#' @export
rastcopula <- function(n, object) {
  if (!methods::is(object, "astcopula")) {
    stop("'object' must be an object of class 'astcopula'.", call. = FALSE)
  }
  ast_sample(n, object@nu)
}

# ---- the base-copula interface ---------------------------------------------

setMethod("basecopula_supported", "astcopula", function(basecopula) TRUE)
setMethod("basecopula_sample", "astcopula", function(basecopula, n) ast_sample(n, basecopula@nu))
setMethod("basecopula_density", "astcopula", function(basecopula, v1, v2) ast_density(v1, v2, basecopula@nu))
setMethod("basecopula_cdf", "astcopula", function(basecopula, v1, v2) ast_cdf(v1, v2, basecopula@nu))
setMethod("basecopula_h", "astcopula", function(basecopula, given, x, given_var) ast_h(given, x, basecopula@nu))
setMethod("basecopula_hinv", "astcopula", function(basecopula, given, q, given_var) ast_hinv(given, q, basecopula@nu))
setMethod("basecopula_fast_hinv", "astcopula", function(basecopula) TRUE)

#' @rdname astcopula
#' @export
setMethod("show", "astcopula", function(object) {
  cat("Absolute spherical t copula, nu = ", format(object@nu), "\n", sep = "")
})

# ---- Kendall's tau and its calibration --------------------------------------
#
# tau = 4 E[C(U, V)] - 1 = 1 - 4 int int h1(u, v) h2(u, v) du dv, a smooth
# double integral over the unit square (Dias, Han and McNeil, equation 34, in
# the form with the h-functions). It is evaluated by tensor tanh-sinh
# quadrature in the complement 1 - u of each variable, with the quantiles
# a = t_nu^-1((1 + u) / 2) computed in the upper tail. The values are
# tau(4) = 0.0994, tau(2) = 0.1894, tau(1) = 1/3 exactly and tau(0.5) = 0.5151,
# and tau(nu) nu -> 4 / pi^2 as nu -> infinity.
ast_tau_exact <- function(nu, N = 400L, tmax = 3.2) {
  t <- seq(-tmax, tmax, length.out = N)
  s <- pi * sinh(t)
  xc <- 1 / (1 + exp(s))
  w <- (t[2] - t[1]) * pi * cosh(t) / ((1 + exp(-s)) * (1 + exp(s)))
  a <- pmin(stats::qt(xc / 2, nu, lower.tail = FALSE), 1e150)
  k <- sqrt((nu + 1) / (nu + a^2))
  h1 <- 1 - 2 * stats::pt(outer(k, a), nu + 1, lower.tail = FALSE)
  h2 <- 1 - 2 * stats::pt(outer(a, k), nu + 1, lower.tail = FALSE)
  1 - 4 * sum(outer(w, w) * h1 * h2)
}

# Monotone cubic interpolation of logit(tau) against log(nu) through the
# tabulated values, with its inverse. The function is almost linear in these
# coordinates (slope -1 at both ends), so the interpolation is very accurate.
ast_tau_interp <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      x <- log(ast_tau_table$nu)
      y <- stats::qlogis(ast_tau_table$tau)
      cache <<- list(
        f = stats::splinefun(x, y, method = "hyman"),
        lo = min(ast_tau_table$nu), hi = max(ast_tau_table$nu),
        tau_lo = max(ast_tau_table$tau), tau_hi = min(ast_tau_table$tau)
      )
    }
    cache
  }
})

#' Kendall's tau of the absolute spherical t copula, and its calibration
#'
#' `astcopula_tau()` gives Kendall's tau of the absolute spherical t copula
#' ([astcopula()]) with `nu` degrees of freedom, and `astcopula_nu()` is its
#' inverse: the degrees of freedom whose Kendall's tau is a given value (the
#' calibration of the copula to a rank correlation). Kendall's tau decreases
#' from 1 to 0 as `nu` increases from 0 to infinity; for example, `tau` is
#' 0.0994, 0.1894, 1/3 and 0.5151 at `nu` = 4, 2, 1 and 0.5.
#'
#' Kendall's tau is a double integral over the unit square of a smooth
#' function of the h-functions, which has no closed form that I know of,
#' except that `tau(1) = 1/3` and that `tau(nu) * nu` tends to `4 / pi^2` as
#' `nu` grows. The default is therefore a table lookup: Kendall's tau was
#' computed once, to about nine digits, at 181 values of `nu` spaced
#' logarithmically between 0.02 and 2000, and the values are interpolated
#' by a monotone cubic spline in the coordinates `log(nu)` and
#' `logit(tau)`, in which the function is nearly linear. Both functions are
#' vectorised and take about a microsecond per value. The error is below
#' `1e-9` in `tau`. Beyond the table, for `nu` above 2000 the asymptotic
#' formula `4 / (pi^2 nu)` (accurate to `1e-14`) is used; `nu` below 0.02,
#' where `tau` exceeds 0.9669, is not supported by the table lookup.
#' `exact = TRUE` evaluates the integral by tanh-sinh quadrature instead
#' (about 10 ms per value), accurate to about `1e-6` for `nu` of 0.05 or
#' more.
#'
#' @param nu degrees of freedom, positive numbers.
#' @param tau values of Kendall's tau in `(0, 1)`.
#' @param exact logical; compute `astcopula_tau()` by quadrature instead of
#'   interpolation?
#' @return A numeric vector.
#' @references
#' Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#' D-vines for stochastic volatility. *Journal of Multivariate Analysis*,
#' **217**, 105695. \doi{10.1016/j.jmva.2026.105695} (the double integral for
#' Kendall's tau is equation 34 of the Supplementary Material).
#' @seealso [astcopula()].
#' @export
#' @examples
#' astcopula_tau(c(0.5, 1, 2, 4, 10))
#' astcopula_nu(c(0.1, 0.3, 0.5))
#' astcopula_tau(astcopula_nu(0.25))
astcopula_tau <- function(nu, exact = FALSE) {
  if (anyNA(nu) || !is.numeric(nu) || any(nu <= 0)) {
    stop("'nu' must be positive numbers.", call. = FALSE)
  }
  if (exact) {
    return(vapply(as.numeric(nu), ast_tau_exact, 0))
  }
  I <- ast_tau_interp()
  if (any(nu < I$lo)) {
    stop(sprintf("the table lookup supports nu >= %g; use exact = TRUE for nu of 0.05 or more.", I$lo),
      call. = FALSE
    )
  }
  out <- 4 / (pi^2 * nu)
  inside <- nu <= I$hi
  out[inside] <- stats::plogis(I$f(log(nu[inside])))
  out
}

#' @rdname astcopula_tau
#' @export
astcopula_nu <- function(tau) {
  if (anyNA(tau) || !is.numeric(tau) || any(tau <= 0 | tau >= 1)) {
    stop("'tau' must be numbers in (0, 1).", call. = FALSE)
  }
  I <- ast_tau_interp()
  if (any(tau > I$tau_lo)) {
    stop(sprintf("the table supports tau <= %.4f (nu >= %g).", I$tau_lo, I$lo), call. = FALSE)
  }
  nu <- 4 / (pi^2 * tau)
  inside <- tau >= I$tau_hi
  if (any(inside)) {
    y <- stats::qlogis(tau[inside])
    # start from the inverse in the table's own coordinates, then polish with
    # Newton steps on the forward interpolant
    x <- stats::approx(I$f(log(ast_tau_table$nu)), log(ast_tau_table$nu), y, rule = 2)$y
    for (i in 1:3) {
      x <- x - (I$f(x) - y) / I$f(x, deriv = 1L)
    }
    nu[inside] <- exp(x)
  }
  nu
}
