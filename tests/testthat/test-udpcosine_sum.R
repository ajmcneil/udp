# Finite orthonormal-cosine-expansion udp transformations: T(u) = F(g(u)) for
# g(u) = sum_j coef[j] * Omega_j(u), Omega_j(u) = sqrt(2)*(-1)^j*cos(j*pi*u).

test_that("chebyshev_coef() matches the known Chebyshev polynomials", {
  known <- list(
    1, c(0, 1), c(-1, 0, 2), c(0, -3, 0, 4),
    c(1, 0, -8, 0, 8), c(0, 5, 0, -20, 0, 16)
  )
  for (d in 0:5) {
    expect_equal(chebyshev_coef(d), known[[d + 1L]])
  }
  # T_j(1) = 1, T_j(-1) = (-1)^j for every degree
  for (d in 0:12) {
    cf <- chebyshev_coef(d)
    expect_equal(polyval(cf, 1), 1)
    expect_equal(polyval(cf, -1), (-1)^d)
  }
})

test_that("udpcosine_sum() validates its arguments", {
  expect_s4_class(udpcosine_sum(c(0, 0, 1)), "udpcosine_sum")
  expect_identical(udpcosine_sum(c(1, -1, 0, 0.5))@degree, 4L)
  expect_error(udpcosine_sum(numeric(0)), "numeric vector")
  expect_error(udpcosine_sum(c(0, 0, 0)), "nonzero entry")
  expect_error(udpcosine_sum(c(1, NA)), "numeric vector")
  expect_error(udpcosine_sum(1, ngrid = 2), "at least 3")
  expect_warning(udpcosine_sum(c(rep(0, 13), 1)), "degree > 12")
})

test_that("trailing zeros are trimmed to find the effective degree", {
  x <- udpcosine_sum(c(0.5, 0, 0.3, 0, 0))
  expect_identical(x@degree, 3L)
  expect_equal(x@coef, c(0.5, 0, 0.3))
})

test_that("a one-hot coefficient reproduces the matching udpcosine, up to 1/sqrt(2)", {
  u <- seq(0, 1, length.out = 1001)
  for (j in 1:6) {
    coef <- rep(0, j)
    coef[j] <- 1 / sqrt(2)
    expect_equal(
      udptrans(udpcosine_sum(coef), u),
      udptrans(udpcosine(j), u),
      tolerance = 1e-4
    )
  }
})

test_that("flipping the sign of every coefficient complements T to 1 - T", {
  u <- seq(0, 1, length.out = 1001)
  coef <- c(0.7, -0.35, 0, 0.5, -0.25)
  a <- udptrans(udpcosine_sum(coef), u)
  b <- udptrans(udpcosine_sum(-coef), u)
  expect_equal(b, 1 - a, tolerance = 1e-8)
})

test_that("udptrans() maps [0, 1] into [0, 1] with no missing values", {
  set.seed(1)
  u <- seq(0, 1, length.out = 1001)
  for (i in 1:8) {
    coef <- round(stats::runif(sample(1:6, 1), -1, 1), 2)
    if (all(coef == 0)) next
    y <- udptrans(udpcosine_sum(coef), u)
    expect_false(anyNA(y))
    expect_true(all(y >= 0 & y <= 1))
  }
})

test_that("udptrans() preserves the shape and attributes of u", {
  u <- ts(seq(0.02, 0.98, length.out = 24), frequency = 4)
  y <- udptrans(udpcosine_sum(c(0.4, 0, -0.6, 0.2)), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("udpcosine_sum transformations preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20260101)
  u <- runif(1e5)
  cases <- list(c(0.8, -0.3, 0.5), c(0, 0, 0, 0.6, -0.4), c(1 / sqrt(2)), c(-0.5, 0.5))
  for (coef in cases) {
    p <- suppressWarnings(
      stats::ks.test(udptrans(udpcosine_sum(coef), u), "punif")$p.value
    )
    expect_gte(p, 0.01)
  }
})

test_that("udpinverse() returns pre-images that invert back through udptrans()", {
  set.seed(2)
  coef <- c(0.6, 0, -0.5, 0.3, 0, -0.2)
  x <- udpcosine_sum(coef)
  v <- runif(200)
  M <- udpinverse(x, v)
  expect_equal(ncol(M), x@degree)
  for (i in seq_along(v)) {
    r <- M[i, ]
    filled <- r[!is.na(r)]
    expect_gte(length(filled), 1L)
    expect_false(is.unsorted(filled))
    expect_true(all(diff(is.na(r)) >= 0)) # NA only trails
    expect_equal(udptrans(x, filled), rep(v[i], length(filled)), tolerance = 5e-3)
  }
})

test_that("udpinverse(prob = TRUE) rows sum to 1 and weight by 1 / |g'|", {
  x <- udpcosine_sum(c(0.5, 0, -0.4, 0.3))
  v <- seq(0.05, 0.95, by = 0.05)
  M <- udpinverse(x, v, prob = TRUE)
  P <- attr(M, "prob")
  expect_identical(dim(P), dim(M))
  expect_identical(is.na(P), is.na(M))
  expect_equal(rowSums(P, na.rm = TRUE), rep(1, length(v)), tolerance = 1e-9)

  present <- !is.na(M)
  w <- 1 / abs(cosine_sum_gderiv(x@cfsD, M[present]))
  expected <- w / ave(w, row(M)[present], FUN = sum)
  expect_equal(P[present], expected, tolerance = 1e-9)
})

test_that("udpinverse() validates v", {
  x <- udpcosine_sum(c(0.5, -0.3, 0.2))
  expect_error(udpinverse(x, 1.2), "in \\[0, 1\\]")
  expect_error(udpinverse(x, c(0.5, NA)), "in \\[0, 1\\]")
  expect_equal(dim(udpinverse(x, numeric(0))), c(0L, x@degree))
})

test_that("udpsi() stochastically inverts udptrans() and recovers a uniform", {
  skip_on_cran()
  set.seed(20260101)
  for (coef in list(c(0.8, -0.3, 0.5), c(0, 0, 0, 0.6, -0.4))) {
    x <- udpcosine_sum(coef)
    u <- runif(1e5)
    back <- udpsi(x, udptrans(x, u), runif(1e5))
    expect_true(all(back >= 0 & back <= 1))
    expect_gte(suppressWarnings(stats::ks.test(back, "punif")$p.value), 0.01)
  }
})

test_that("udpderiv() matches a central finite difference away from special points", {
  set.seed(4)
  coef <- c(0.6, 0, -0.5, 0.3, 0, -0.2)
  x <- udpcosine_sum(coef)
  tp_x <- chebyshev_turnpoints(x@cfsD)
  tp_u <- if (length(tp_x)) sort(acos(pmin(pmax(tp_x, -1), 1)) / pi) else numeric(0)
  h <- 1e-6
  u <- runif(300)
  u <- u[vapply(u, function(u0) all(abs(u0 - tp_u) > 0.01) && u0 > 0.01 && u0 < 0.99, logical(1))]
  u <- u[1:80]
  d <- udpderiv(x, u)
  fd <- (udptrans(x, u + h) - udptrans(x, u - h)) / (2 * h)
  expect_equal(d, fd, tolerance = 1e-2)
})

test_that("udpderiv() at u = 0 / u = 1 matches a one-sided finite difference", {
  h <- 1e-6
  for (coef in list(c(0.7, -0.35, 0, 0.5, -0.25), c(0.4, 0.3), c(-0.6, 0.2, 0.1, -0.3))) {
    x <- udpcosine_sum(coef)
    fd0 <- as.numeric((udptrans(x, h) - udptrans(x, 0)) / h)
    fd1 <- as.numeric((udptrans(x, 1) - udptrans(x, 1 - h)) / h)
    expect_equal(udpderiv(x, 0), fd0, tolerance = 1e-2)
    expect_equal(udpderiv(x, 1), fd1, tolerance = 1e-2)
  }
})

test_that("udpderiv() boundary multiplicity handles a pure single-term reduction", {
  # coef = c(0, 0, 1) reduces to a pure cosine term (like udpcosine(3)): every
  # trough/peak of the periodic g shares its extreme value with u = 0 or
  # u = 1, so the naive multiplicity-1 boundary formula is wrong by a factor
  # of 3 here -- this is the case that originally caught the bug.
  x <- udpcosine_sum(c(0, 0, 1))
  h <- 1e-6
  fd0 <- as.numeric((udptrans(x, h) - udptrans(x, 0)) / h)
  fd1 <- as.numeric((udptrans(x, 1) - udptrans(x, 1 - h)) / h)
  expect_equal(udpderiv(x, 0), fd0, tolerance = 1e-2)
  expect_equal(udpderiv(x, 1), fd1, tolerance = 1e-2)
  expect_equal(abs(udpderiv(x, 0)), 3, tolerance = 1e-2)
  expect_equal(abs(udpderiv(x, 1)), 3, tolerance = 1e-2)
})

test_that("pcoincide() lies in (0, 1] and agrees with the default integrator", {
  x <- udpcosine_sum(c(0.5, 0, -0.4, 0.3))
  p <- pcoincide(x)
  expect_gt(p, 0)
  expect_lte(p, 1)
  expect_equal(p, integrate_collision(x), tolerance = 1e-4)
})

test_that("plot() runs, including for a degree with no turning points", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  expect_no_error(plot(udpcosine_sum(c(1)))) # degree 1: tp is empty
  expect_no_error(plot(udpcosine_sum(c(0.7, -0.35, 0, 0.5, -0.25))))
  expect_no_error(plot(udpcosine_sum(c(0, 0.6, 0, -0.4)), embellish = "bw"))
  expect_no_error(plot(udpcosine_sum(c(0.3, -0.6)), n = 200, embellish = "none"))
})

test_that("plot() draws horizontal lines at T(0) and T(1) when non-trivial", {
  # u = 0/1 are always A-partition members by definition (every Omega_j has
  # zero derivative at both), but their T-values need not be trivial or
  # match any interior turning point -- a regression check for the gap where
  # filtering 0/1 out of the vertical lines also, incorrectly, filtered them
  # out of the horizontal ones.
  x <- udpcosine_sum(c(0.7, -0.35, 0, 0.5, -0.25))
  t0 <- as.numeric(udptrans(x, 0))
  t1 <- as.numeric(udptrans(x, 1))
  expect_true(t0 > 1e-6 && t0 < 1 - 1e-6) # sanity: genuinely non-trivial
  expect_true(t1 > 1e-6 && t1 < 1 - 1e-6)

  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  hlines <- numeric(0)
  local_mocked_bindings(
    segments = function(x0, y0, x1, y1, ...) {
      if (identical(x0, 0) && identical(x1, 1)) hlines <<- c(hlines, y0)
      invisible(NULL)
    },
    .package = "udp"
  )
  plot(x, embellish = "colour")
  expect_true(any(abs(hlines - t0) < 1e-6))
  expect_true(any(abs(hlines - t1) < 1e-6))
})
