test_that("both branches of vsi() invert the v-transform: V(vsi(x, v)) == v", {
  v <- seq(0.05, 0.95, length.out = 19)
  n <- length(v)
  for (case in vtransform_list()) {
    down <- vsi(case$x, v, Z = rep(0, n))
    up   <- vsi(case$x, v, Z = rep(1, n))
    expect_equal(vtrans(case$x, down), v, tolerance = 1e-6)
    expect_equal(vtrans(case$x, up), v, tolerance = 1e-6)
  }
})

test_that("Z = 0 selects the lower pre-image, Z = 1 the upper one", {
  v <- seq(0.05, 0.95, length.out = 19)
  n <- length(v)
  for (case in vtransform_list()) {
    lower <- vinverse(case$x, v)
    expect_equal(vsi(case$x, v, Z = rep(0, n)), lower, tolerance = 1e-8)
    expect_equal(vsi(case$x, v, Z = rep(1, n)), v + lower, tolerance = 1e-8)
  }
})

test_that("vsi() output lies in [0, 1] for a random Z", {
  set.seed(42)
  v <- runif(200)
  for (case in vtransform_list()) {
    out <- vsi(case$x, v)
    expect_true(all(out >= 0 & out <= 1))
    expect_length(out, length(v))
  }
})

test_that("vsi() errors when Z and v have different lengths", {
  expect_error(
    vsi(vsymmetric(), c(0.2, 0.4, 0.6), Z = c(0, 1)),
    "same length"
  )
})

test_that("vsi() is reproducible under a fixed seed and preserves ts attributes", {
  v <- ts(runif(50), frequency = 12)

  set.seed(1); a <- vsi(v2p(delta = 0.4, kappa = 1.2), v)
  set.seed(1); b <- vsi(v2p(delta = 0.4, kappa = 1.2), v)
  expect_identical(a, b)
  expect_s3_class(a, "ts")
  expect_identical(tsp(a), tsp(v))
})
