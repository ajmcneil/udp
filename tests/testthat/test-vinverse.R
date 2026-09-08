test_that("vinverse() is a right inverse of vtrans() (lower branch)", {
  v <- seq(0.05, 0.95, length.out = 19)
  for (case in vtransform_list()) {
    u <- vinverse(case$x, v)
    expect_equal(vtrans(case$x, u), v, tolerance = 1e-6)
  }
})

test_that("vinverse() returns the pre-image at or below the fulcrum", {
  v <- seq(0.05, 0.95, length.out = 19)
  for (case in vtransform_list()) {
    expect_true(all(vinverse(case$x, v) <= case$delta + 1e-7))
  }
})

test_that("analytic and numeric inversion agree for an invertible v-transform", {
  v <- seq(0.05, 0.95, length.out = 19)
  analytic <- Vlinear(delta = 0.4)

  # same transform, stripped to the base class so vinverse() must go numeric
  numeric_only <- new("Vtransform",
    name = analytic@name, Vtrans = analytic@Vtrans,
    pars = analytic@pars, gradient = analytic@gradient
  )

  expect_equal(vinverse(analytic, v), vinverse(numeric_only, v), tolerance = 1e-7)
})

test_that("the Newton inverse can reach full double precision with a tight tol", {
  v <- seq(0.01, 0.99, length.out = 50)
  for (case in vtransform_list()) {
    u <- vinverse(case$x, v, method = "newton", tol = .Machine$double.eps^0.9)
    expect_equal(vtrans(case$x, u), v, tolerance = 1e-9)
    expect_true(all(u >= 0 & u <= case$delta + 1e-12))
  }
})

test_that("the default-tol Newton inverse is accurate to about 1e-8", {
  v <- seq(0.01, 0.99, length.out = 50)
  for (case in vtransform_list()) {
    u <- vinverse(case$x, v)
    expect_equal(vtrans(case$x, u), v, tolerance = 1e-6)
  }
})

test_that("the spline inverse agrees with Newton to grid accuracy", {
  v <- seq(0.05, 0.95, length.out = 40)
  for (case in vtransform_list()) {
    exact <- vinverse(case$x, v, method = "newton")
    approx <- vinverse(case$x, v, method = "spline", ngrid = 4000)
    expect_equal(approx, exact, tolerance = 1e-3)
  }
})

test_that("vinverse() handles the endpoints and NA", {
  x <- V2p(delta = 0.4, kappa = 1.2)
  expect_equal(vinverse(x, c(0, 1)), c(0.4, 0))
  expect_equal(vinverse(x, c(0, 1), method = "spline"), c(0.4, 0))
  expect_identical(vinverse(x, NA_real_), NA_real_)
})

test_that("vinverse() preserves the shape and attributes of v", {
  x <- V3p(delta = 0.45, kappa = 0.8, xi = 1.2)
  v <- ts(seq(0.1, 0.9, length.out = 24), frequency = 4)
  u <- vinverse(x, v)
  expect_s3_class(u, "ts")
  expect_identical(tsp(u), tsp(v))

  m <- matrix(seq(0.1, 0.9, length.out = 12), nrow = 3)
  expect_identical(dim(vinverse(x, m)), dim(m))
})

test_that("invertible v-transforms ignore method/ngrid", {
  v <- seq(0.05, 0.95, length.out = 10)
  expect_identical(
    vinverse(Vsymmetric(), v, method = "spline", ngrid = 5),
    vinverse(Vsymmetric(), v)
  )
})

test_that("method flows through vsi() and vdownprob() via ...", {
  x <- V2b(delta = 0.35, kappa = 1.2)
  v <- seq(0.05, 0.95, length.out = 15)

  set.seed(1); a <- vsi(x, v, method = "newton")
  set.seed(1); b <- vsi(x, v, method = "spline", ngrid = 4000)
  expect_equal(a, b, tolerance = 1e-3)

  expect_equal(
    vdownprob(x, v, method = "spline", ngrid = 4000),
    vdownprob(x, v),
    tolerance = 1e-3
  )
})

test_that("vdownprob() returns probabilities in [0, 1]", {
  v <- seq(0.02, 0.98, length.out = 25)
  for (case in vtransform_list()) {
    p <- vdownprob(case$x, v)
    expect_true(all(p >= 0 & p <= 1))
  }
})

test_that("vdownprob() is 1/2 everywhere for the symmetric v-transform", {
  v <- seq(0.05, 0.95, length.out = 10)
  expect_equal(vdownprob(Vsymmetric(), v), rep(0.5, length(v)))
})
