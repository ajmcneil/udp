test_that("pcoincide() is exactly 1/2 for the symmetric v-transform", {
  expect_identical(pcoincide(vsymmetric()), 0.5)
})

test_that("pcoincide() uses the closed form for the linear v-transform", {
  for (d in c(0.2, 0.4, 0.5, 0.75)) {
    expect_equal(pcoincide(vlinear(delta = d)), d^2 + (1 - d)^2)
  }

  # the closed form agrees with numerical integration of the definition
  d <- 0.35
  x <- vlinear(delta = d)
  numeric <- d^2 + (1 - d)^2 +
    2 * stats::integrate(function(v) (vdownprob(x, v) - d)^2, 0, 1)$value
  expect_equal(pcoincide(x), numeric, tolerance = 1e-8)
})

test_that("pcoincide() matches a direct Monte Carlo estimate", {
  skip_on_cran()
  set.seed(2024)
  x <- v2p(delta = 0.4, kappa = 1.3)

  u <- runif(2e5)
  v <- udptrans(x, u)
  u_back <- udpsi(x, v)
  # a "coincidence" is landing back on the same pre-image; the numeric
  # inverse is accurate to ~1e-8, distinct pre-images differ by O(0.1)
  mc <- mean(abs(u_back - u) < 1e-6)

  expect_equal(pcoincide(x), mc, tolerance = 0.02)
})

test_that("pcoincide() is exactly 1 for a shuffle", {
  expect_identical(pcoincide(shuffle(c(2, 1, 3))), 1)
  expect_identical(pcoincide(shuffle(c(3, 1, 2), signs = c(1, -1, 1))), 1)
})

test_that("pcoincide() is exactly 1 / degree for a cosine udp transformation", {
  for (k in 1:8) {
    expect_identical(pcoincide(udpcosine(k)), 1 / k)
  }
})

test_that("pcoincide() returns a single unnamed probability in (0, 1]", {
  cases <- c(
    lapply(vtransform_list(), `[[`, "x"),
    list(shuffle(c(2, 1, 3)), udpcosine(4), udplegendre(3), udplegendre(6))
  )
  for (x in cases) {
    p <- pcoincide(x)
    expect_length(p, 1)
    expect_null(names(p))
    expect_gt(p, 0)
    expect_lte(p, 1 + 1e-9)
  }
})

test_that("the default udp method integrates to the closed forms", {
  general <- selectMethod("pcoincide", "udp")
  expect_equal(general(udpcosine(3)), 1 / 3, tolerance = 1e-6)
  expect_equal(general(udpcosine(5)), 1 / 5, tolerance = 1e-6)
  expect_equal(general(shuffle(c(3, 1, 2))), 1, tolerance = 1e-6)

  x <- v2p(delta = 0.4, kappa = 1.3)
  expect_equal(general(x), pcoincide(x), tolerance = 1e-4)
})

test_that("pcoincide() is 1 for degree 1 and 1/2 for degree 2 Legendre", {
  expect_equal(pcoincide(udplegendre(1)), 1, tolerance = 1e-8)
  expect_equal(pcoincide(udplegendre(2)), 0.5, tolerance = 1e-6)
})

test_that("pcoincide() for udplegendre matches a Monte Carlo estimate", {
  skip_on_cran()
  set.seed(2025)
  for (j in c(3, 4, 6, 9)) {
    x <- udplegendre(j)
    # a "coincidence" is landing back on the same monotone branch of L_j.
    # Compare branch indices rather than |u_back - u|: a distance threshold
    # counts wrong pre-images that happen to sit close to u near the turning
    # points of L_j, a bias that grows with the degree.
    turns <- sort(Re(polyroot(x@cfsD)))
    turns <- turns[turns > 0 & turns < 1]
    branch <- function(z) findInterval(z, turns)
    u <- runif(3e5)
    u_back <- udpsi(x, udptrans(x, u))
    mc <- mean(branch(u_back) == branch(u))
    # tolerance covers Monte Carlo noise plus pcoincide()'s own ~1e-3 accuracy
    # for udplegendre (the spline-interpolated F_j)
    expect_equal(pcoincide(x), mc, tolerance = 0.02)
  }
})
