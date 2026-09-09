# Shifted-Legendre udp transformations: T(u) = F_j(L_j(u)).

# exact F_j(y) for testing, straight from the definition
leg_cdf_exact <- function(degree, y) {
  cfs <- slegendre_coef(degree)
  lb <- legendre_lbound(cfs, poly_deriv_coef(cfs), degree)
  vapply(y, legendre_measure, numeric(1), coef = cfs, lbound = lb)
}

test_that("slegendre_coef() matches the known shifted Legendre polynomials", {
  known <- list(
    1, c(-1, 2), c(1, -6, 6), c(-1, 12, -30, 20),
    c(1, -20, 90, -140, 70), c(-1, 30, -210, 560, -630, 252)
  )
  for (d in 0:5) {
    expect_equal(slegendre_coef(d), known[[d + 1L]])
  }
  # endpoint identities L_j(1) = 1, L_j(0) = (-1)^j
  for (d in 0:12) {
    cf <- slegendre_coef(d)
    expect_equal(polyval(cf, 1), 1)
    expect_equal(polyval(cf, 0), (-1)^d)
  }
})

test_that("poly_deriv_coef() and polyval() agree with a finite difference", {
  cf <- slegendre_coef(6)
  cfD <- poly_deriv_coef(cf)
  x <- c(0.05, 0.3, 0.55, 0.8, 0.95)
  fd <- (polyval(cf, x + 1e-6) - polyval(cf, x - 1e-6)) / 2e-6
  expect_equal(polyval(cfD, x), fd, tolerance = 1e-6)
})

test_that("udplegendre() validates its arguments", {
  expect_s4_class(udplegendre(3), "udplegendre")
  expect_identical(udplegendre(4)@degree, 4L)
  expect_error(udplegendre(2.5), "positive integer")
  expect_error(udplegendre(0), "positive integer")
  expect_error(udplegendre(c(2, 3)), "positive integer")
  expect_error(udplegendre(NA), "positive integer")
  expect_error(udplegendre(5, ngrid = 2), "at least 3")
  expect_warning(udplegendre(14), "degree > 12")
})

test_that("degree 1 is the identity and degree 2 is the tent map", {
  u <- seq(0, 1, length.out = 101)
  expect_equal(udptrans(udplegendre(1), u), u)
  expect_equal(udptrans(udplegendre(2), u), abs(2 * u - 1))
})

test_that("udptrans() maps [0, 1] into [0, 1] and matches exact F_j(L_j(u))", {
  u <- seq(0, 1, length.out = 1001)
  for (d in 1:9) {
    x <- udplegendre(d)
    y <- udptrans(x, u)
    expect_false(anyNA(y))
    expect_true(all(y >= 0 & y <= 1))
    expect_equal(y, leg_cdf_exact(d, polyval(x@cfs, u)), tolerance = 2e-3)
  }
})

test_that("udptrans() preserves the shape and attributes of u", {
  u <- ts(seq(0.02, 0.98, length.out = 24), frequency = 4)
  y <- udptrans(udplegendre(5), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("shifted-Legendre transformations preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (d in 1:8) {
    p <- suppressWarnings(stats::ks.test(udptrans(udplegendre(d), u), "punif")$p.value)
    expect_gte(p, 0.01)
  }
})

test_that("udpinverse() returns a degree-column matrix of pre-images", {
  set.seed(1)
  for (d in c(1, 3, 4, 6, 8)) {
    x <- udplegendre(d)
    v <- runif(200)
    M <- udpinverse(x, v)
    expect_equal(dim(M), c(200L, d))
    for (i in seq_along(v)) {
      r <- M[i, ]
      filled <- r[!is.na(r)]
      expect_gte(length(filled), 1L)
      expect_false(is.unsorted(filled))
      expect_true(all(diff(is.na(r)) >= 0)) # NA only trails
      expect_equal(udptrans(x, filled), rep(v[i], length(filled)), tolerance = 5e-3)
    }
  }
})

test_that("udpinverse() pre-image count varies with v", {
  # a degree-5 map does not always have 5 real pre-images
  x <- udplegendre(5)
  set.seed(2)
  counts <- rowSums(!is.na(udpinverse(x, runif(500))))
  expect_true(any(counts < 5L))
  expect_true(all(counts >= 1L & counts <= 5L))
})

test_that("udpinverse(prob = TRUE) weights pre-images by 1 / |L_j'|", {
  x <- udplegendre(4)
  v <- seq(0.05, 0.95, by = 0.05)
  M <- udpinverse(x, v, prob = TRUE)
  P <- attr(M, "prob")
  expect_identical(dim(P), dim(M))
  expect_identical(is.na(P), is.na(M))
  expect_equal(rowSums(P, na.rm = TRUE), rep(1, length(v)), tolerance = 1e-9)

  present <- !is.na(M)
  w <- 1 / abs(polyval(poly_deriv_coef(x@cfs), M[present]))
  expected <- w / ave(w, row(M)[present], FUN = sum)
  expect_equal(P[present], expected, tolerance = 1e-9)
})

test_that("udpinverse() validates v", {
  x <- udplegendre(3)
  expect_error(udpinverse(x, 1.2), "in \\[0, 1\\]")
  expect_error(udpinverse(x, c(0.5, NA)), "in \\[0, 1\\]")
  expect_equal(dim(udpinverse(x, numeric(0))), c(0L, 3L))
})

test_that("udpsi() stochastically inverts udptrans() and recovers a uniform", {
  skip_on_cran()
  set.seed(20240908)
  for (d in 1:8) {
    x <- udplegendre(d)
    u <- runif(1e5)
    back <- udpsi(x, udptrans(x, u), runif(1e5))
    expect_true(all(back >= 0 & back <= 1))
    expect_gte(suppressWarnings(stats::ks.test(back, "punif")$p.value), 0.01)
  }
})

test_that("udpsi() on a shifted-Legendre map is reproducible and keeps attributes", {
  x <- udplegendre(4)
  v <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  set.seed(3)
  a <- udpsi(x, v)
  set.seed(3)
  b <- udpsi(x, v)
  expect_identical(a, b)
  expect_s3_class(a, "ts")
  expect_identical(tsp(a), tsp(v))
})

test_that("plot() runs for a udplegendre", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  for (d in 1:7) {
    expect_no_error(plot(udplegendre(d)))
  }
  expect_no_error(plot(udplegendre(5), n = 200, main = "degree 5", ylab = "y"))
})
