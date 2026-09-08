test_that("invertible constructors return vtransformi objects", {
  for (ctor in list(vsymmetric(), vlinear(delta = 0.4))) {
    expect_s4_class(ctor, "vtransformi")
    expect_s4_class(ctor, "vtransform")
  }
})

test_that("parametric constructors return (non-invertible) vtransform objects", {
  for (ctor in list(
    v2p(delta = 0.4, kappa = 1.2),
    v2b(delta = 0.4, kappa = 1.2),
    v3p(delta = 0.4, kappa = 1.2, xi = 1.1),
    v3b(delta = 0.4, kappa = 1.2, xi = 1.1)
  )) {
    expect_s4_class(ctor, "vtransform")
    expect_false(is(ctor, "vtransformi"))
  }
})

test_that("constructors record their own name", {
  expect_identical(vsymmetric()@name, "vsymmetric")
  expect_identical(vlinear()@name, "vlinear")
  expect_identical(v2p()@name, "v2p")
  expect_identical(v2b()@name, "v2b")
  expect_identical(v3p()@name, "v3p")
  expect_identical(v3b()@name, "v3b")
})

test_that("constructors reject a fulcrum outside (0, 1)", {
  for (ctor in list(vlinear, v2p, v2b, v3p, v3b)) {
    expect_error(ctor(delta = 0), "in \\(0, 1\\)")
    expect_error(ctor(delta = 1), "in \\(0, 1\\)")
    expect_error(ctor(delta = -0.1), "in \\(0, 1\\)")
    expect_error(ctor(delta = c(0.3, 0.4)), "single number")
  }
})

test_that("coef() returns the named parameter vector", {
  expect_equal(coef(v2p(delta = 0.4, kappa = 1.2)), c(delta = 0.4, kappa = 1.2))
  expect_equal(
    coef(v3p(delta = 0.3, kappa = 0.8, xi = 1.1)),
    c(delta = 0.3, kappa = 0.8, xi = 1.1)
  )
})

test_that("show() method prints the name and parameters", {
  expect_output(show(v2p(delta = 0.4, kappa = 1.2)), "v2p")
  expect_output(show(v2p(delta = 0.4, kappa = 1.2)), "delta")
})
