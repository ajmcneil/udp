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
    expect_true(all(vinverse(case$x, v) <= case$delta + 1e-8))
  }
})

test_that("analytic and numeric inversion agree for an invertible v-transform", {
  v <- seq(0.05, 0.95, length.out = 19)
  analytic <- Vlinear(delta = 0.4)

  # same transform, but stripped to the base class so vinverse() must use uniroot
  numeric_only <- new("Vtransform",
    name = analytic@name, Vtrans = analytic@Vtrans,
    pars = analytic@pars, gradient = analytic@gradient
  )

  expect_equal(vinverse(analytic, v), vinverse(numeric_only, v), tolerance = 1e-6)
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
