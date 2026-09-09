test_that("udpcosine() validates degree and stores it as an integer", {
  x <- udpcosine(3)
  expect_s4_class(x, "udpcosine")
  expect_identical(x@degree, 3L)

  expect_error(udpcosine(2.5), "positive integer")
  expect_error(udpcosine(0), "positive integer")
  expect_error(udpcosine(-1), "positive integer")
  expect_error(udpcosine(c(1, 2)), "positive integer")
  expect_error(udpcosine(NA), "positive integer")
})

test_that("udptrans() matches known closed forms", {
  u <- seq(0, 1, length.out = 101)
  expect_equal(udptrans(udpcosine(1), u), u)            # degree 1 is the identity
  expect_equal(udptrans(udpcosine(2), u), abs(2 * u - 1)) # degree 2 is the tent
})

test_that("udptrans() maps [0, 1] into [0, 1]", {
  u <- seq(0, 1, length.out = 2001)
  for (d in 1:8) {
    y <- udptrans(udpcosine(d), u)
    expect_false(anyNA(y))
    expect_true(all(y >= 0 & y <= 1))
  }
})

test_that("udptrans() is piecewise linear with slope of magnitude degree", {
  for (d in 1:7) {
    x <- udpcosine(d)
    k <- seq_len(d)
    u1 <- (k - 0.75) / d # two points strictly inside piece k = [(k-1)/d, k/d]
    u2 <- (k - 0.25) / d
    slope <- (udptrans(x, u2) - udptrans(x, u1)) / (u2 - u1)
    expect_equal(abs(slope), rep(d, d), tolerance = 1e-6)
  }
})

test_that("udptrans() preserves the shape and attributes of u", {
  u <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  y <- udptrans(udpcosine(3), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("cosine udp transformations preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (d in 1:6) {
    p <- suppressWarnings(stats::ks.test(udptrans(udpcosine(d), u), "punif")$p.value)
    expect_gte(p, 0.01)
  }
})

test_that("udpcosinverse() returns a list aligned with v", {
  x <- udpcosine(4)
  res <- udpcosinverse(x, c(0.2, 0.6))
  expect_type(res, "list")
  expect_length(res, 2L)

  named <- udpcosinverse(x, c(a = 0.3, b = 0.7))
  expect_named(named, c("a", "b"))
})

test_that("udpcosinverse() gives 'degree' roots for v in (0, 1), each a pre-image", {
  set.seed(1)
  for (d in 1:8) {
    x <- udpcosine(d)
    v <- runif(30)
    res <- udpcosinverse(x, v)
    for (j in seq_along(v)) {
      r <- res[[j]]
      expect_length(r, d)
      expect_false(is.unsorted(r))
      expect_true(all(r >= 0 & r <= 1))
      expect_equal(udptrans(x, r), rep(v[j], d), tolerance = 1e-9)
    }
  }
})

test_that("udpcosinverse() collapses shared roots at v = 0 and v = 1", {
  for (d in 1:6) {
    x <- udpcosine(d)
    r0 <- udpcosinverse(x, 0)[[1]]
    r1 <- udpcosinverse(x, 1)[[1]]

    expect_length(r0, ceiling(d / 2))
    expect_length(r1, floor(d / 2) + 1)
    expect_equal(udptrans(x, r0), rep(0, length(r0)), tolerance = 1e-9)
    expect_equal(udptrans(x, r1), rep(1, length(r1)), tolerance = 1e-9)
    expect_false(anyDuplicated(r0) > 0)
    expect_false(anyDuplicated(r1) > 0)
  }
})

test_that("udpcosinverse() matches a known small case", {
  expect_equal(
    udpcosinverse(udpcosine(3), c(0, 0.4, 1)),
    list(c(0, 2 / 3), c(2 / 15, 8 / 15, 4 / 5), c(1 / 3, 1))
  )
  expect_equal(udpcosinverse(udpcosine(1), c(0, 0.25, 1)), list(0, 0.25, 1))
})

test_that("udpcosinverse() validates v", {
  x <- udpcosine(3)
  expect_error(udpcosinverse(x, 1.5), "in \\[0, 1\\]")
  expect_error(udpcosinverse(x, -0.01), "in \\[0, 1\\]")
  expect_error(udpcosinverse(x, c(0.5, NA)), "in \\[0, 1\\]")
})

test_that("udpsi() returns one pre-image per v in (0, 1)", {
  set.seed(1)
  for (d in 1:6) {
    x <- udpcosine(d)
    v <- runif(1000)
    Z <- runif(1000)
    u <- udpsi(x, v, Z)

    expect_length(u, 1000L)
    expect_true(all(u >= 0 & u <= 1))
    expect_equal(udptrans(x, u), v, tolerance = 1e-9)
  }
})

test_that("udpsi() uses Z as a categorical quantile over the roots", {
  set.seed(2)
  x <- udpcosine(5)
  v <- runif(500)
  roots <- udpcosinverse(x, v)

  expect_equal(
    udpsi(x, v, rep(1e-9, 500)),
    vapply(roots, `[`, numeric(1), 1)
  )
  expect_equal(
    udpsi(x, v, rep(1, 500)),
    vapply(roots, `[`, numeric(1), 5)
  )
  expect_equal(
    udpsi(x, v, rep(0.5, 500)), # ceiling(0.5 * 5) = 3
    vapply(roots, `[`, numeric(1), 3)
  )
})

test_that("udpsi() at v = 0 and v = 1 samples among the shared roots", {
  set.seed(1)
  for (d in 1:6) {
    x <- udpcosine(d)
    v <- rep(c(0, 1), each = 40)
    u <- udpsi(x, v, runif(length(v)))
    expect_true(all(u >= 0 & u <= 1))
    expect_equal(udptrans(x, u), v, tolerance = 1e-9)
    expect_true(all(u[v == 0] %in% udpcosinverse(x, 0)[[1]]))
    expect_true(all(u[v == 1] %in% udpcosinverse(x, 1)[[1]]))
  }
})

test_that("udpsi() validates its arguments and edge inputs", {
  x <- udpcosine(3)
  expect_error(udpsi(x, c(0.3, 0.4), Z = 0.5), "same length")
  expect_error(udpsi(x, 1.2, 0.5), "in \\[0, 1\\]")
  expect_error(udpsi(x, c(0.5, NA), c(0.1, 0.1)), "in \\[0, 1\\]")
  expect_length(udpsi(x, numeric(0), numeric(0)), 0L)
})

test_that("udpsi() is reproducible and preserves attributes", {
  x <- udpcosine(3)
  set.seed(7)
  a <- udpsi(x, seq(0.1, 0.9, by = 0.1))
  set.seed(7)
  b <- udpsi(x, seq(0.1, 0.9, by = 0.1))
  expect_identical(a, b)

  v <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  out <- udpsi(x, v)
  expect_s3_class(out, "ts")
  expect_identical(tsp(out), tsp(v))
})

test_that("the stochastic inverse recovers a uniform input", {
  skip_on_cran()
  set.seed(20240908)
  for (d in 1:6) {
    x <- udpcosine(d)
    u <- runif(1e5)
    back <- udpsi(x, udptrans(x, u), runif(1e5))
    expect_gte(suppressWarnings(stats::ks.test(back, "punif")$p.value), 0.01)
  }
})

test_that("plot() runs for a udpcosine", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  for (d in 1:5) {
    expect_no_error(plot(udpcosine(d)))
  }
  expect_no_error(plot(udpcosine(3), main = "degree 3", ylab = "y"))
})
