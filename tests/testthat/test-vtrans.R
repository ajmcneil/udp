test_that("v-transforms fix the defining points V(0)=1, V(delta)=0, V(1)=1", {
  for (case in vtransform_list()) {
    expect_equal(udptrans(case$x, 0), 1, tolerance = 1e-8)
    expect_equal(udptrans(case$x, case$delta), 0, tolerance = 1e-8)
    expect_equal(udptrans(case$x, 1), 1, tolerance = 1e-8)
  }
})

test_that("v-transforms map [0, 1] into [0, 1]", {
  u <- c(0, u_grid(), 1)
  for (case in vtransform_list()) {
    vv <- udptrans(case$x, u)
    expect_true(all(vv >= 0 & vv <= 1))
    expect_length(vv, length(u))
  }
})

test_that("v-transforms are decreasing then increasing about the fulcrum", {
  for (case in vtransform_list()) {
    below <- seq(0.01, case$delta - 0.01, length.out = 15)
    above <- seq(case$delta + 0.01, 0.99, length.out = 15)
    expect_true(all(diff(udptrans(case$x, below)) < 0))
    expect_true(all(diff(udptrans(case$x, above)) > 0))
  }
})

test_that("udptrans() is vectorised and preserves length", {
  x <- v2p(delta = 0.4, kappa = 1.2)
  expect_length(udptrans(x, u_grid(30)), 30)
})

test_that("vsymmetric and vlinear match their closed forms", {
  u <- u_grid()
  expect_equal(udptrans(vsymmetric(), u), abs(2 * u - 1))

  delta <- 0.4
  expected <- abs(u / delta - 1) * ((delta / (1 - delta))^(u > delta))
  expect_equal(udptrans(vlinear(delta = delta), u), expected)
})
