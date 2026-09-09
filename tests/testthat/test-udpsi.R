test_that("both branches of udpsi() invert the v-transform: V(udpsi(x, v)) == v", {
  v <- seq(0.05, 0.95, length.out = 19)
  n <- length(v)
  for (case in vtransform_list()) {
    down <- udpsi(case$x, v, Z = rep(0, n))
    up   <- udpsi(case$x, v, Z = rep(1, n))
    expect_equal(udptrans(case$x, down), v, tolerance = 1e-6)
    expect_equal(udptrans(case$x, up), v, tolerance = 1e-6)
  }
})

test_that("Z = 0 selects the lower pre-image, Z = 1 the upper one", {
  v <- seq(0.05, 0.95, length.out = 19)
  n <- length(v)
  for (case in vtransform_list()) {
    lower <- vinverse(case$x, v)
    expect_equal(udpsi(case$x, v, Z = rep(0, n)), lower, tolerance = 1e-8)
    expect_equal(udpsi(case$x, v, Z = rep(1, n)), v + lower, tolerance = 1e-8)
  }
})

test_that("udpsi() output lies in [0, 1] for a random Z", {
  set.seed(42)
  v <- runif(200)
  for (case in vtransform_list()) {
    out <- udpsi(case$x, v)
    expect_true(all(out >= 0 & out <= 1))
    expect_length(out, length(v))
  }
})

test_that("udpsi() errors when Z and v have different lengths", {
  expect_error(
    udpsi(vsymmetric(), c(0.2, 0.4, 0.6), Z = c(0, 1)),
    "same length"
  )
})

test_that("udpsi() is reproducible under a fixed seed and preserves ts attributes", {
  v <- ts(runif(50), frequency = 12)

  set.seed(1); a <- udpsi(v2p(delta = 0.4, kappa = 1.2), v)
  set.seed(1); b <- udpsi(v2p(delta = 0.4, kappa = 1.2), v)
  expect_identical(a, b)
  expect_s3_class(a, "ts")
  expect_identical(tsp(a), tsp(v))
})

test_that("udpsi() on a shuffle is the deterministic inverse and ignores Z", {
  s <- shuffle(c(3, 1, 2), signs = c(1, -1, 1))
  v <- c(0.1, 0.5, 0.9)
  expect_identical(udpsi(s, v), shinverse(s, v))
  expect_identical(udpsi(s, v, Z = rep(0, 3)), shinverse(s, v))
})

test_that("udpsi() on a shuffle does not consume the RNG stream", {
  s <- shuffle(c(2, 1, 3))
  v <- c(0.2, 0.4, 0.7)
  set.seed(1); r1 <- runif(1)
  set.seed(1); invisible(udpsi(s, v)); r2 <- runif(1)
  expect_identical(r1, r2)
})

test_that("udpsi() inverts a cosine udp transformation", {
  x <- udpcosine(3)
  set.seed(1)
  u <- runif(200)
  back <- udpsi(x, udptrans(x, u))
  expect_true(all(back >= 0 & back <= 1))
  expect_equal(udptrans(x, back), udptrans(x, u), tolerance = 1e-9)
})
