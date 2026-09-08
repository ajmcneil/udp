test_that("v-transforms fix the defining points V(0)=1, V(delta)=0, V(1)=1", {
  for (case in vtransform_list()) {
    expect_equal(vtrans(case$x, 0), 1, tolerance = 1e-8)
    expect_equal(vtrans(case$x, case$delta), 0, tolerance = 1e-8)
    expect_equal(vtrans(case$x, 1), 1, tolerance = 1e-8)
  }
})

test_that("v-transforms map [0, 1] into [0, 1]", {
  u <- c(0, u_grid(), 1)
  for (case in vtransform_list()) {
    vv <- vtrans(case$x, u)
    expect_true(all(vv >= 0 & vv <= 1))
    expect_length(vv, length(u))
  }
})

test_that("v-transforms are decreasing then increasing about the fulcrum", {
  for (case in vtransform_list()) {
    below <- seq(0.01, case$delta - 0.01, length.out = 15)
    above <- seq(case$delta + 0.01, 0.99, length.out = 15)
    expect_true(all(diff(vtrans(case$x, below)) < 0))
    expect_true(all(diff(vtrans(case$x, above)) > 0))
  }
})

test_that("vtrans() is vectorised and preserves length", {
  x <- V2p(delta = 0.4, kappa = 1.2)
  expect_length(vtrans(x, u_grid(30)), 30)
})

test_that("Vsymmetric and Vlinear match their closed forms", {
  u <- u_grid()
  expect_equal(vtrans(Vsymmetric(), u), abs(2 * u - 1))

  delta <- 0.4
  expected <- abs(u / delta - 1) * ((delta / (1 - delta))^(u > delta))
  expect_equal(vtrans(Vlinear(delta = delta), u), expected)
})

test_that("Vdegenerate is the identity", {
  u <- u_grid()
  expect_equal(vtrans(Vdegenerate(), u), u)
  expect_equal(vinverse(Vdegenerate(), u), u)
})

test_that("delta = 0 and delta = 1 collapse the parametric v-transforms", {
  u <- c(0.2, 0.5, 0.8)
  for (ctor in list(V2p, V2b, function(d) V3p(d), function(d) V3b(d))) {
    expect_equal(vtrans(ctor(0), u), u)
    expect_equal(vtrans(ctor(1), u), 1 - u)
  }
})
