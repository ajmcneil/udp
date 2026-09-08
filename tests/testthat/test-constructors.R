test_that("invertible constructors return VtransformI objects", {
  for (ctor in list(Vsymmetric(), Vdegenerate(), Vlinear(delta = 0.4))) {
    expect_s4_class(ctor, "VtransformI")
    expect_s4_class(ctor, "Vtransform")
  }
})

test_that("parametric constructors return (non-invertible) Vtransform objects", {
  for (ctor in list(
    V2p(delta = 0.4, kappa = 1.2),
    V2b(delta = 0.4, kappa = 1.2),
    V3p(delta = 0.4, kappa = 1.2, xi = 1.1),
    V3b(delta = 0.4, kappa = 1.2, xi = 1.1)
  )) {
    expect_s4_class(ctor, "Vtransform")
    expect_false(is(ctor, "VtransformI"))
  }
})

test_that("constructors record their own name", {
  expect_identical(Vsymmetric()@name, "Vsymmetric")
  expect_identical(Vdegenerate()@name, "Vdegenerate")
  expect_identical(Vlinear()@name, "Vlinear")
  expect_identical(V2p()@name, "V2p")
  expect_identical(V2b()@name, "V2b")
  expect_identical(V3p()@name, "V3p")
  expect_identical(V3b()@name, "V3b")
})

test_that("coef() returns the named parameter vector", {
  expect_equal(coef(V2p(delta = 0.4, kappa = 1.2)), c(delta = 0.4, kappa = 1.2))
  expect_equal(
    coef(V3p(delta = 0.3, kappa = 0.8, xi = 1.1)),
    c(delta = 0.3, kappa = 0.8, xi = 1.1)
  )
})

test_that("show() method prints the name and parameters", {
  expect_output(show(V2p(delta = 0.4, kappa = 1.2)), "V2p")
  expect_output(show(V2p(delta = 0.4, kappa = 1.2)), "delta")
})
