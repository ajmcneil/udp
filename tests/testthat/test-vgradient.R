test_that("vgradient() agrees with a finite-difference derivative of vtrans()", {
  h <- 1e-6
  for (case in vtransform_list()) {
    below <- seq(0.05, case$delta - 0.05, length.out = 8)
    above <- seq(case$delta + 0.05, 0.95, length.out = 8)
    u <- c(below, above)
    fd <- (vtrans(case$x, u + h) - vtrans(case$x, u - h)) / (2 * h)
    expect_equal(vgradient(case$x, u), fd, tolerance = 1e-4)
  }
})

test_that("the gradient is negative below the fulcrum and positive above", {
  for (case in vtransform_list()) {
    below <- seq(0.05, case$delta - 0.05, length.out = 10)
    above <- seq(case$delta + 0.05, 0.95, length.out = 10)
    expect_true(all(vgradient(case$x, below) < 0))
    expect_true(all(vgradient(case$x, above) > 0))
  }
})

test_that("Vsymmetric and Vlinear gradients match their closed forms", {
  expect_equal(vgradient(Vsymmetric(), c(0.2, 0.8)), c(-2, 2))

  delta <- 0.4
  expect_equal(
    vgradient(Vlinear(delta = delta), c(0.2, 0.8)),
    c(-1 / delta, 1 / (1 - delta))
  )
})
