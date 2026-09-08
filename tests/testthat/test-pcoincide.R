test_that("pcoincide() is exactly 1/2 for the symmetric v-transform", {
  expect_identical(pcoincide(Vsymmetric()), 0.5)
})

test_that("pcoincide() returns a single unnamed probability in (0, 1]", {
  for (case in vtransform_list()) {
    p <- pcoincide(case$x)
    expect_length(p, 1)
    expect_null(names(p))
    expect_gt(p, 0)
    expect_lte(p, 1)
  }
})

test_that("pcoincide() matches a direct Monte Carlo estimate", {
  skip_on_cran()
  set.seed(2024)
  x <- V2p(delta = 0.4, kappa = 1.3)

  u <- runif(2e5)
  v <- vtrans(x, u)
  u_back <- vsi(x, v)
  # a "coincidence" is landing back on the same pre-image; the numeric
  # inverse is accurate to ~1e-8, distinct pre-images differ by O(0.1)
  mc <- mean(abs(u_back - u) < 1e-6)

  expect_equal(pcoincide(x), mc, tolerance = 0.02)
})
