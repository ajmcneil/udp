# Finite orthonormal-Legendre-expansion udp transformations: T(u) = F(g(u))
# for g(u) = sum_j coef[j] * P_j(u), P_j(u) = sqrt(2j + 1) * L_j(u).

test_that("udplegendre_sum() validates its arguments", {
  expect_s4_class(udplegendre_sum(c(0, 0, 1)), "udplegendre_sum")
  expect_identical(udplegendre_sum(c(1, -1, 0, 0.5))@degree, 4L)
  expect_error(udplegendre_sum(numeric(0)), "numeric vector")
  expect_error(udplegendre_sum(c(0, 0, 0)), "nonzero entry")
  expect_error(udplegendre_sum(c(1, NA)), "numeric vector")
  expect_error(udplegendre_sum(1, ngrid = 2), "at least 3")
  expect_warning(udplegendre_sum(c(rep(0, 13), 1)), "degree > 12")
})

test_that("trailing zeros are trimmed to find the effective degree", {
  x <- udplegendre_sum(c(0.5, 0, 0.3, 0, 0))
  expect_identical(x@degree, 3L)
  expect_equal(x@coef, c(0.5, 0, 0.3))
})

test_that("a one-hot positive coefficient reproduces the matching udplegendre", {
  u <- seq(0, 1, length.out = 1001)
  for (j in 1:6) {
    coef <- rep(0, j)
    coef[j] <- 1
    expect_equal(
      udptrans(udplegendre_sum(coef), u),
      udptrans(udplegendre(j), u),
      tolerance = 1e-4
    )
  }
})

test_that("flipping the sign of every coefficient complements T to 1 - T", {
  u <- seq(0, 1, length.out = 1001)
  coef <- c(0.7, -0.35, 0, 0.5, -0.25)
  a <- udptrans(udplegendre_sum(coef), u)
  b <- udptrans(udplegendre_sum(-coef), u)
  expect_equal(b, 1 - a, tolerance = 1e-8)
})

test_that("udptrans() maps [0, 1] into [0, 1] with no missing values", {
  set.seed(1)
  u <- seq(0, 1, length.out = 1001)
  for (i in 1:8) {
    coef <- round(stats::runif(sample(1:6, 1), -1, 1), 2)
    if (all(coef == 0)) next
    y <- udptrans(udplegendre_sum(coef), u)
    expect_false(anyNA(y))
    expect_true(all(y >= 0 & y <= 1))
  }
})

test_that("udptrans() preserves the shape and attributes of u", {
  u <- ts(seq(0.02, 0.98, length.out = 24), frequency = 4)
  y <- udptrans(udplegendre_sum(c(0.4, 0, -0.6, 0.2)), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("udplegendre_sum transformations preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20260101)
  u <- runif(1e5)
  cases <- list(c(0.8, -0.3, 0.5), c(0, 0, 0, 0.6, -0.4), c(1), c(-0.5, 0.5))
  for (coef in cases) {
    p <- suppressWarnings(
      stats::ks.test(udptrans(udplegendre_sum(coef), u), "punif")$p.value
    )
    expect_gte(p, 0.01)
  }
})

test_that("udpinverse() returns pre-images that invert back through udptrans()", {
  set.seed(2)
  coef <- c(0.6, 0, -0.5, 0.3, 0, -0.2)
  x <- udplegendre_sum(coef)
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

test_that("udpinverse(prob = TRUE) weights pre-images by 1 / |g'| and rows sum to 1", {
  x <- udplegendre_sum(c(0.5, 0, -0.4, 0.3))
  v <- seq(0.05, 0.95, by = 0.05)
  M <- udpinverse(x, v, prob = TRUE)
  P <- attr(M, "prob")
  expect_identical(dim(P), dim(M))
  expect_identical(is.na(P), is.na(M))
  expect_equal(rowSums(P, na.rm = TRUE), rep(1, length(v)), tolerance = 1e-9)

  present <- !is.na(M)
  w <- 1 / abs(polyval(x@cfsD, M[present]))
  expected <- w / ave(w, row(M)[present], FUN = sum)
  expect_equal(P[present], expected, tolerance = 1e-9)
})

test_that("udpinverse() validates v", {
  x <- udplegendre_sum(c(0.5, -0.3, 0.2))
  expect_error(udpinverse(x, 1.2), "in \\[0, 1\\]")
  expect_error(udpinverse(x, c(0.5, NA)), "in \\[0, 1\\]")
  expect_equal(dim(udpinverse(x, numeric(0))), c(0L, x@degree))
})

test_that("udpsi() stochastically inverts udptrans() and recovers a uniform", {
  skip_on_cran()
  set.seed(20260101)
  for (coef in list(c(0.8, -0.3, 0.5), c(0, 0, 0, 0.6, -0.4))) {
    x <- udplegendre_sum(coef)
    u <- runif(1e5)
    back <- udpsi(x, udptrans(x, u), runif(1e5))
    expect_true(all(back >= 0 & back <= 1))
    expect_gte(suppressWarnings(stats::ks.test(back, "punif")$p.value), 0.01)
  }
})

test_that("pcoincide() lies in (0, 1] and agrees with the default integrator", {
  x <- udplegendre_sum(c(0.5, 0, -0.4, 0.3))
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
  expect_no_error(plot(udplegendre_sum(c(1)))) # degree 1: tp is empty
  expect_no_error(plot(udplegendre_sum(c(0.7, -0.35, 0, 0.5, -0.25))))
  expect_no_error(plot(udplegendre_sum(c(0, 0.6, 0, -0.4)), embellish = "bw"))
  expect_no_error(plot(udplegendre_sum(c(0.3, -0.6)), n = 200, embellish = "none"))
})

test_that("plot() draws a horizontal line at T(1) when it is a non-trivial T-partition value", {
  # coef = c(0.3, -0.2): u = 1 is not a turning point and does not share a
  # turning point's critical value, but u = 0/1 are always A-partition
  # members by definition, so T(1) must still appear as a T-partition
  # horizontal line -- a regression check for the gap where filtering 0/1 out
  # of the vertical lines also, incorrectly, filtered them out of the
  # horizontal ones.
  x <- udplegendre_sum(c(0.3, -0.2))
  target <- as.numeric(udptrans(x, 1))
  expect_true(target > 1e-6 && target < 1 - 1e-6) # sanity: genuinely non-trivial

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
  expect_true(any(abs(hlines - target) < 1e-6))
})
