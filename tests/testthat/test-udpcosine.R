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

test_that("udpcostrans() matches known closed forms", {
  u <- seq(0, 1, length.out = 101)
  expect_equal(udpcostrans(udpcosine(1), u), u)            # degree 1 is the identity
  expect_equal(udpcostrans(udpcosine(2), u), abs(2 * u - 1)) # degree 2 is the tent
})

test_that("udpcostrans() maps [0, 1] into [0, 1]", {
  u <- seq(0, 1, length.out = 2001)
  for (d in 1:8) {
    y <- udpcostrans(udpcosine(d), u)
    expect_false(anyNA(y))
    expect_true(all(y >= 0 & y <= 1))
  }
})

test_that("udpcostrans() is piecewise linear with slope of magnitude degree", {
  for (d in 1:7) {
    x <- udpcosine(d)
    k <- seq_len(d)
    u1 <- (k - 0.75) / d # two points strictly inside piece k = [(k-1)/d, k/d]
    u2 <- (k - 0.25) / d
    slope <- (udpcostrans(x, u2) - udpcostrans(x, u1)) / (u2 - u1)
    expect_equal(abs(slope), rep(d, d), tolerance = 1e-6)
  }
})

test_that("udpcostrans() preserves the shape and attributes of u", {
  u <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  y <- udpcostrans(udpcosine(3), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("cosine udp transformations preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (d in 1:6) {
    p <- suppressWarnings(stats::ks.test(udpcostrans(udpcosine(d), u), "punif")$p.value)
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
      expect_equal(udpcostrans(x, r), rep(v[j], d), tolerance = 1e-9)
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
    expect_equal(udpcostrans(x, r0), rep(0, length(r0)), tolerance = 1e-9)
    expect_equal(udpcostrans(x, r1), rep(1, length(r1)), tolerance = 1e-9)
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
