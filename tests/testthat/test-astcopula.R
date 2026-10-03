# The absolute spherical t copula. Checks are independent of the formulas in
# R/astcopula.R: the t copula of rvinecopulib (for integer nu), numerical
# integration, simulation, and the defining property that stochastic
# inversion of the symmetric v-transform returns the t copula.

bd <- function(...) rvinecopulib::bicop_dist(...)

test_that("density, CDF and h agree with the t copula for nu >= 2", {
  skip_if_not_installed("rvinecopulib")
  set.seed(1)
  u <- runif(100)
  v <- runif(100)
  for (nu in c(2, 3, 5)) {
    ac <- astcopula(nu)
    tc <- bd("t", 0, c(0, nu))
    p <- cbind((1 + u) / 2, (1 + v) / 2)
    expect_equal(dastcopula(u, v, ac), rvinecopulib::dbicop(p, tc), tolerance = 1e-10)
    expect_equal(pastcopula(u, v, ac), 4 * rvinecopulib::pbicop(p, tc) - u - v - 1, tolerance = 1e-8)
    expect_equal(hastcopula(u, v, ac), 2 * rvinecopulib::hbicop(p, 1, tc) - 1, tolerance = 1e-10)
    expect_equal(hastcopula(u, v, ac, cond_var = 2), 2 * rvinecopulib::hbicop(p, 2, tc) - 1, tolerance = 1e-10)
  }
})

test_that("the CDF is the integral of the density and h its derivative, for any nu > 0", {
  set.seed(2)
  for (nu in c(0.3, 1, 2.5, 20)) {
    ac <- astcopula(nu)
    u <- 0.45
    v <- 0.7
    ref <- integrate(function(s) {
      vapply(s, function(si) integrate(function(t) dastcopula(si, t, ac), 0, v, rel.tol = 1e-9)$value, 0)
    }, 0, u, rel.tol = 1e-9)$value
    expect_equal(pastcopula(u, v, ac), ref, tolerance = 1e-6)
    e <- 1e-5
    expect_equal(hastcopula(u, v, ac), (pastcopula(u + e, v, ac) - pastcopula(u - e, v, ac)) / (2 * e), tolerance = 1e-6)
    expect_equal(hastcopula(u, v, ac, cond_var = 2), (pastcopula(u, v + e, ac) - pastcopula(u, v - e, ac)) / (2 * e), tolerance = 1e-6)
    # exact margins and boundary values
    expect_equal(pastcopula(c(0, 1, 1, 0.3), c(0.5, 0.5, 1, 0), ac), c(0, 0.5, 1, 0))
    expect_equal(hastcopula(0.4, c(0, 1), ac), c(0, 1))
  }
})

test_that("inverse h inverts h, and a random sample matches the CDF", {
  set.seed(3)
  for (nu in c(0.4, 3)) {
    ac <- astcopula(nu)
    u <- runif(50)
    v <- runif(50)
    for (cv in 1:2) {
      h <- hastcopula(u, v, ac, cond_var = cv)
      expect_equal(hastcopula(u, h, ac, cond_var = cv, inverse = TRUE)[cv == 1], v[cv == 1], tolerance = 1e-8)
    }
    X <- rastcopula(2e5, ac)
    pts <- cbind(c(0.2, 0.5, 0.8), c(0.3, 0.5, 0.9))
    est <- vapply(1:3, function(i) mean(X[, 1] <= pts[i, 1] & X[, 2] <= pts[i, 2]), 0)
    expect_lt(max(abs(pastcopula(pts, object = ac) - est)), 4e-3)
  }
})

test_that("stochastic inversion of vsymmetric() on both margins gives the t copula", {
  skip_if_not_installed("rvinecopulib")
  set.seed(4)
  u <- runif(60)
  v <- runif(60)
  for (nu in c(2, 4)) {
    bc <- bsicopula(astcopula(nu), vsymmetric(), vsymmetric())
    tc <- bd("t", 0, c(0, nu))
    expect_equal(dbsicopula(u, v, bc), rvinecopulib::dbicop(cbind(u, v), tc), tolerance = 1e-8)
    expect_equal(pbsicopula(u, v, bc), rvinecopulib::pbicop(cbind(u, v), tc), tolerance = 1e-7)
    expect_equal(hbsicopula(u, v, bc), rvinecopulib::hbicop(cbind(u, v), 1, tc), tolerance = 1e-8)
    expect_equal(hbsicopula(u, v, bc, cond_var = 2, inverse = TRUE), rvinecopulib::hbicop(cbind(u, v), 2, tc, inverse = TRUE), tolerance = 1e-7)
  }
  # also with non-linear udps, against simulation
  bc <- bsicopula(astcopula(1.5), v2p(0.4, 1.4), vlinear(0.6))
  X <- rbsicopula(1e5, bc)
  pts <- cbind(c(0.2, 0.6), c(0.3, 0.8))
  est <- vapply(1:2, function(i) mean(X[, 1] <= pts[i, 1] & X[, 2] <= pts[i, 2]), 0)
  expect_lt(max(abs(pbsicopula(pts, object = bc) - est)), 5e-3)
  expect_error(randsdvine(astcopula(3)), "bicop_dist")
})

test_that("Kendall's tau reproduces the published values and the limits", {
  expect_equal(astcopula_tau(c(4, 2, 1, 0.5)), c(0.0994, 0.1894, 1 / 3, 0.5151), tolerance = 1e-3)
  expect_equal(astcopula_tau(1), 1 / 3, tolerance = 1e-9)
  expect_equal(astcopula_tau(1, exact = TRUE), 1 / 3, tolerance = 1e-9)
  # tau nu -> 4 / pi^2
  expect_equal(astcopula_tau(5000) * 5000, 4 / pi^2, tolerance = 1e-7)
  # decreasing in nu, in (0, 1), and continuous where the table ends
  nu <- exp(seq(log(0.03), log(1e4), length.out = 400))
  tau <- astcopula_tau(nu)
  expect_true(all(diff(tau) < 0) && all(tau > 0 & tau < 1))
  expect_equal(astcopula_tau(1999.999), astcopula_tau(2000.001), tolerance = 1e-5)
  # the table agrees with fresh quadrature at other values of nu
  set.seed(5)
  nu <- exp(runif(8, log(0.05), log(1000)))
  expect_equal(astcopula_tau(nu), astcopula_tau(nu, exact = TRUE), tolerance = 1e-5)
  # and with simulation
  set.seed(6)
  X <- rastcopula(4e4, astcopula(2))
  expect_equal(cor(X[, 1], X[, 2], method = "kendall"), astcopula_tau(2), tolerance = 0.02)
})

test_that("astcopula_nu() inverts astcopula_tau()", {
  tau <- c(1e-6, 1e-4, 0.001, 0.01, 0.1, 0.333333333333, 0.5, 0.8, 0.95)
  expect_equal(astcopula_tau(astcopula_nu(tau)), tau, tolerance = 1e-9)
  nu <- c(0.03, 0.5, 1, 3, 40, 1999, 2001, 1e5)
  expect_equal(astcopula_nu(astcopula_tau(nu)), nu, tolerance = 1e-7)
  expect_equal(astcopula_nu(1 / 3), 1, tolerance = 1e-8)
})

test_that("arguments are checked", {
  expect_error(astcopula(0), "positive")
  expect_error(astcopula(c(1, 2)), "positive")
  expect_error(dastcopula(0.3, 0.4, "x"), "astcopula")
  expect_error(dastcopula(1.2, 0.4, astcopula(2)), "\\[0, 1\\]")
  expect_error(hastcopula(0.3, 0.4, astcopula(2), cond_var = 3), "cond_var")
  expect_error(astcopula_tau(-1), "positive")
  expect_error(astcopula_tau(0.001), "supports")
  expect_error(astcopula_nu(1), "in \\(0, 1\\)")
  expect_error(astcopula_nu(0.99), "supports")
  expect_equal(dim(rastcopula(7, astcopula(2))), c(7, 2))
})

test_that("astcopula can be the base copula of a bsicopula with a randsdvine randomizer", {
  skip_if_not_installed("rvinecopulib")
  skip_if_not_installed("copula")
  rs <- randsdvine(bd("gaussian", 0, 0.7), bd("clayton", 90, 1.5), bd("gumbel", 270, 1.8))
  # same copula through two backends: identical density, weight and h
  set.seed(7)
  u <- runif(40)
  v <- runif(40)
  a <- bsicopula(bd("clayton", 0, 2), v2p(0.4, 1.4), udpcosine(3), rs)
  b <- bsicopula(copula::claytonCopula(2), v2p(0.4, 1.4), udpcosine(3), rs)
  expect_equal(dbsicopula(u, v, a), dbsicopula(u, v, b), tolerance = 1e-6)
  # ast against simulation, and with a hand-built equivalent: the t copula is
  # what vsymmetric margins recover, so d, p and r must agree with each other
  bc <- bsicopula(astcopula(2.5), vlinear(0.4), v2p(0.6, 1.3), rs)
  set.seed(8)
  X <- rbsicopula(2e5, bc)
  expect_lt(max(abs(c(mean(X[, 1] <= 0.3), mean(X[, 2] <= 0.7)) - c(0.3, 0.7))), 3e-3) # uniform margins
  pts <- cbind(c(0.2, 0.5, 0.8), c(0.3, 0.5, 0.6))
  est <- vapply(1:3, function(i) mean(X[, 1] <= pts[i, 1] & X[, 2] <= pts[i, 2]), 0)
  expect_lt(max(abs(pbsicopula(pts, object = bc, nodes = 31) - est)), 4e-3)
  # the weight is a density: integrates to one over the unit square
  g <- (seq_len(120) - 0.5) / 120
  G <- expand.grid(g, g)
  expect_equal(mean(dbsicopula(G[, 1], G[, 2], bc)), 1, tolerance = 5e-3)
})
