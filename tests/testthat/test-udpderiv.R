# udpderiv(): the derivative of udptrans(), including its convention (the
# left derivative, or the right derivative at u = 0) at non-smooth points.

test_that("udpderiv() matches a finite-difference derivative away from kinks", {
  h <- 1e-6
  cases <- list(
    list(x = shuffle(c(3, 1, 2), signs = c(1, -1, 1)),
         u = c(0.05, 0.25, 0.5, 0.61, 0.95), tol = 1e-6),
    list(x = udpcosine(4), u = c(0.03, 0.2, 0.4, 0.55, 0.9), tol = 1e-6),
    list(x = vsymmetric(), u = c(0.1, 0.3, 0.7, 0.9), tol = 1e-5),
    list(x = v2p(delta = 0.4, kappa = 1.3), u = c(0.05, 0.2, 0.6, 0.95), tol = 1e-5),
    list(x = udplegendre(4), u = c(0.05, 0.15, 0.45, 0.85, 0.95), tol = 1e-3)
  )
  for (case in cases) {
    fd <- (udptrans(case$x, case$u + h) - udptrans(case$x, case$u - h)) / (2 * h)
    expect_equal(udpderiv(case$x, case$u), fd, tolerance = case$tol)
  }
})

test_that("udpderiv() preserves the shape and attributes of u", {
  u <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  for (x in list(shuffle(c(2, 1, 3)), udpcosine(3), vsymmetric(), udplegendre(4))) {
    d <- udpderiv(x, u)
    expect_s3_class(d, "ts")
    expect_identical(tsp(d), tsp(u))
  }
})

test_that("shuffle: udpderiv() is the strip sign, left strip at a boundary", {
  s <- shuffle(c(3, 1, 2), signs = c(1, -1, 1))
  # strips: [0, 1/3] sign 1, [1/3, 2/3] sign -1, [2/3, 1] sign 1
  expect_equal(udpderiv(s, c(0, 1 / 3, 2 / 3, 1)), c(1, 1, -1, 1))
})

test_that("udpcosine: udpderiv() is +-degree, left piece at a kink", {
  x <- udpcosine(3)
  # pieces: (0,1/3) increasing, (1/3,2/3) decreasing, (2/3,1) increasing
  expect_equal(udpderiv(x, c(0, 1 / 3, 2 / 3, 1)), c(3, 3, -3, 3))
})

test_that("vtransform: udpderiv() is vgradient()", {
  for (case in vtransform_list()) {
    expect_identical(udpderiv(case$x, u_grid()), vgradient(case$x, u_grid()))
  }
  expect_equal(udpderiv(vlinear(delta = 0.4), 0.4), -1 / 0.4)
})

test_that("udplegendre: degree 1 and 2 derivatives match their closed forms", {
  expect_equal(udpderiv(udplegendre(1), c(0, 0.3, 0.7, 1)), rep(1, 4))
  # T(u) = |2u - 1|; the corner at u = 0.5 takes the left derivative -2
  expect_equal(
    udpderiv(udplegendre(2), c(0, 0.3, 0.5, 0.7, 1)),
    c(-2, -2, -2, 2, 2)
  )
})

test_that("udplegendre: the general formula agrees with the turning-point special case", {
  # points just outside the turning-point catchment go through the general
  # f_j(L_j(u)) * L_j'(u) formula; they should already sit close to the exact
  # corner value returned at the turning point itself
  eps <- 2e-4
  for (deg in 3:8) {
    x <- udplegendre(deg)
    tp <- legendre_turnpoints(x@cfsD)
    exact <- udpderiv(x, tp)
    left <- udpderiv(x, pmax(tp - eps, 0))
    right <- udpderiv(x, pmin(tp + eps, 1))
    expect_equal(left, exact, tolerance = 1e-2)
    expect_equal(right, -exact, tolerance = 1e-2)
  }
})

test_that("udplegendre: turning-point slope magnitude accounts for mirror pairs", {
  # even degree is symmetric about u = 1/2 (L_j(1-u) = L_j(u)): the centre
  # turning point is unpaired (magnitude 2), every other one is paired with
  # its mirror and shares a critical value (magnitude 4)
  for (deg in c(4, 6, 8)) {
    x <- udplegendre(deg)
    tp <- legendre_turnpoints(x@cfsD)
    mag <- abs(udpderiv(x, tp))
    expect_equal(mag[abs(tp - 0.5) < 1e-6], 2)
    expect_true(all(mag[abs(tp - 0.5) >= 1e-6] == 4))
  }
  # odd degree is antisymmetric (L_j(1-u) = -L_j(u)): no coincident critical
  # values, every turning point keeps magnitude 2
  for (deg in c(3, 5, 7)) {
    x <- udplegendre(deg)
    tp <- legendre_turnpoints(x@cfsD)
    expect_equal(abs(udpderiv(x, tp)), rep(2, length(tp)))
  }
})

test_that("udplegendre: udpderiv() blows up at a transversal crossing of a turning-point value", {
  # degree 3: the turning point at u ~ 0.2113 (a local max of L_3) and its
  # transversal re-crossing at u ~ 0.9273 share T(0.2113); T' should blow up
  # approaching the crossing from the singular side
  x <- udplegendre(3)
  tp <- legendre_turnpoints(x@cfsD)[1]
  y <- polyval(x@cfs, tp)
  u1 <- legendre_realroots(x@cfs, y)
  u1 <- u1[abs(u1 - tp) > 1e-3]
  expect_length(u1, 1L)
  near <- udpderiv(x, u1 - 10^-(3:6))
  expect_true(all(diff(abs(near)) > 0)) # grows without bound approaching u1
  expect_gt(abs(near[4]), 100)
})
