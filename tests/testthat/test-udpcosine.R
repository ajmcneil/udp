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
