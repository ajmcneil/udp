# pbsicopula() and hbsicopula(): the CDF and h-functions of a bsicopula with
# independent randomizers. Needs rvinecopulib throughout.
#
# The checks are independent of the implementation: the paper's closed forms
# for linear v-transforms (Proposition S3), transcribed as printed; the
# defining relations (h is the derivative of the CDF, the CDF is the integral
# of the density, the CDF matches simulation); and the base copulas include
# 90 and 270 degree rotations, which are not exchangeable, so a swapped
# argument cannot go unnoticed.

bd <- function(...) rvinecopulib::bicop_dist(...)

# Proposition S3 of the supplement to Dias, Han and McNeil, as printed. V is
# the linear v-transform; s(u, d) = d^I{u <= d} (d - 1)^I{u > d}; Cstar,
# h1star = dCstar/du and h2star = dCstar/dv are the base copula's functions.
s3 <- local({
  V <- function(u, d) udptrans(vlinear(d), u)
  s <- function(u, d) ifelse(u <= d, d, d - 1)
  list(
    C = function(u, v, d1, d2, cop) {
      s(u, d1) * s(v, d2) * rvinecopulib::pbicop(cbind(V(u, d1), V(v, d2)), cop) + d1 * v + d2 * u - d1 * d2
    },
    h1 = function(u, v, d1, d2, cop) {
      -s(v, d2) * rvinecopulib::hbicop(cbind(V(u, d1), V(v, d2)), 1, cop) + d2
    },
    h2 = function(u, v, d1, d2, cop) {
      -s(u, d1) * rvinecopulib::hbicop(cbind(V(u, d1), V(v, d2)), 2, cop) + d1
    },
    h1inv = function(u, v, d1, d2, cop) {
      -s(v, d2) * rvinecopulib::hbicop(cbind(V(u, d1), V(v, d2)), 1, cop, inverse = TRUE) + d2
    },
    h2inv = function(u, v, d1, d2, cop) {
      -s(u, d1) * rvinecopulib::hbicop(cbind(V(u, d1), V(v, d2)), 2, cop, inverse = TRUE) + d1
    }
  )
})

# bsicopula models with piecewise-linear margins of every kind, over base
# copulas that are not exchangeable.
linear_models <- function() {
  T3 <- udpzigzag(widths = c(3, 4, 3))
  list(
    `clayton90 | vlinear x zigzag` = bsicopula(bd("clayton", 90, 2), vlinear(0.4), T3),
    `gumbel270 | zigzag(up) x cosine(3)` = bsicopula(bd("gumbel", 270, 2), udpzigzag(widths = c(1, 2, 1), up = TRUE), udpcosine(3)),
    `bb1 90 | vsymmetric x cosine(4)` = bsicopula(bd("bb1", 90, c(0.7, 1.5)), vsymmetric(), udpcosine(4)),
    `clayton270 | shuffle x zigzag` = bsicopula(bd("clayton", 270, 2), shuffle(c(3, 1, 2), signs = c(1, -1, 1)), T3),
    `gaussian | udpid x udpflip` = bsicopula(bd("gaussian", 0, 0.6), udpid(), udpflip()),
    `frank | shuffle x shuffle` = bsicopula(bd("frank", 0, 5), shuffle(c(2, 3, 1), signs = c(-1, 1, 1)), shuffle(c(4, 1, 3, 2), signs = c(1, 1, -1, -1)))
  )
}

test_that("pbsicopula() and hbsicopula() reproduce Proposition S3 for linear v-transforms", {
  skip_if_not_installed("rvinecopulib")
  set.seed(1)
  u <- runif(400)
  v <- runif(400)
  for (cop in list(bd("clayton", 90, 2), bd("gumbel", 270, 2), bd("gaussian", 0, 0.5), bd("joe", 0, 2))) {
    # each of the four combinations of (u <= d1, v <= d2) occurs
    d1 <- 0.35
    d2 <- 0.7
    o <- bsicopula(cop, vlinear(d1), vlinear(d2))
    expect_equal(pbsicopula(u, v, o), s3$C(u, v, d1, d2, cop), tolerance = 1e-9)
    expect_equal(hbsicopula(u, v, o), s3$h1(u, v, d1, d2, cop), tolerance = 1e-9)
    expect_equal(hbsicopula(u, v, o, cond_var = 2), s3$h2(u, v, d1, d2, cop), tolerance = 1e-9)
    expect_equal(hbsicopula(u, v, o, inverse = TRUE), s3$h1inv(u, v, d1, d2, cop), tolerance = 1e-9)
    expect_equal(hbsicopula(u, v, o, cond_var = 2, inverse = TRUE), s3$h2inv(u, v, d1, d2, cop), tolerance = 1e-9)
  }
})

test_that("the CDF has uniform margins and the h-functions are monotone, bounded and invertible", {
  skip_if_not_installed("rvinecopulib")
  set.seed(2)
  u <- runif(300)
  g <- seq(0, 1, length.out = 201)
  for (nm in names(linear_models())) {
    o <- linear_models()[[nm]]
    # C(u, 0) = 0 = C(0, v), C(u, 1) = u, C(1, v) = v
    expect_equal(pbsicopula(u, 0, o), rep(0, 300), info = nm)
    expect_equal(pbsicopula(0, u, o), rep(0, 300), info = nm)
    expect_equal(pbsicopula(u, 1, o), u, tolerance = 1e-9, info = nm)
    expect_equal(pbsicopula(1, u, o), u, tolerance = 1e-9, info = nm)
    for (cv in 1:2) {
      # h runs from 0 to 1 and increases in its target argument
      h <- if (cv == 1) hbsicopula(0.37, g, o, cond_var = 1) else hbsicopula(g, 0.37, o, cond_var = 2)
      expect_equal(h[c(1, 201)], c(0, 1), info = nm)
      expect_true(all(diff(h) >= -1e-9), info = nm)
      # h(inverse(p)) = p
      p <- runif(300)
      uinv <- if (cv == 1) hbsicopula(u, p, o, cond_var = 1, inverse = TRUE) else hbsicopula(p, u, o, cond_var = 2, inverse = TRUE)
      back <- if (cv == 1) hbsicopula(u, uinv, o, cond_var = 1) else hbsicopula(uinv, u, o, cond_var = 2)
      expect_equal(back, p, tolerance = 1e-7, info = paste(nm, "cond_var", cv))
    }
  }
})

test_that("inverse(h(u)) recovers u when the base copula has positive density", {
  skip_if_not_installed("rvinecopulib")
  set.seed(3)
  u <- runif(300)
  w <- runif(300)
  for (cop in list(bd("gaussian", 0, 0.6), bd("frank", 0, 5), bd("gaussian", 0, -0.7))) {
    for (o in list(
      bsicopula(cop, vlinear(0.4), udpzigzag(widths = c(3, 4, 3))),
      bsicopula(cop, udpcosine(3), shuffle(c(3, 1, 2), signs = c(1, -1, 1)))
    )) {
      h1 <- hbsicopula(u, w, o)
      expect_equal(hbsicopula(u, h1, o, inverse = TRUE), w, tolerance = 1e-6)
      h2 <- hbsicopula(w, u, o, cond_var = 2)
      expect_equal(hbsicopula(h2, u, o, cond_var = 2, inverse = TRUE), w, tolerance = 1e-6)
    }
  }
})

test_that("h is the partial derivative of the CDF, for both conditioning variables", {
  skip_if_not_installed("rvinecopulib")
  set.seed(4)
  u <- runif(200, 0.02, 0.98)
  v <- runif(200, 0.02, 0.98)
  d <- 1e-6
  for (nm in names(linear_models())) {
    o <- linear_models()[[nm]]
    dC1 <- (pbsicopula(u + d, v, o) - pbsicopula(u - d, v, o)) / (2 * d)
    dC2 <- (pbsicopula(u, v + d, o) - pbsicopula(u, v - d, o)) / (2 * d)
    expect_equal(hbsicopula(u, v, o, cond_var = 1), dC1, tolerance = 1e-5, info = nm)
    expect_equal(hbsicopula(u, v, o, cond_var = 2), dC2, tolerance = 1e-5, info = nm)
  }
})

test_that("the CDF matches simulation and the integral of the density", {
  skip_if_not_installed("rvinecopulib")
  set.seed(5)
  q <- rbind(c(0.2, 0.3), c(0.5, 0.5), c(0.8, 0.4), c(0.35, 0.9))
  for (nm in c("clayton90 | vlinear x zigzag", "clayton270 | shuffle x zigzag", "bb1 90 | vsymmetric x cosine(4)")) {
    o <- linear_models()[[nm]]
    U <- rbsicopula(2e5, o)
    sim <- apply(q, 1, function(z) mean(U[, 1] <= z[1] & U[, 2] <= z[2]))
    expect_lt(max(abs(pbsicopula(q[, 1], q[, 2], o) - sim)), 0.004)
  }
  # integral of the density over [0, u1] x [0, u2], by the midpoint rule
  o <- linear_models()[["clayton90 | vlinear x zigzag"]]
  g1 <- (seq_len(300) - 0.5) / 300 * 0.5
  g2 <- (seq_len(300) - 0.5) / 300 * 0.5
  G <- expand.grid(g1, g2)
  expect_equal(pbsicopula(0.5, 0.5, o), mean(dbsicopula(G[, 1], G[, 2], o)) * 0.25, tolerance = 5e-3)
})

test_that("h integrates the density when the conditioning margin is not linear", {
  skip_if_not_installed("rvinecopulib")
  # only the margin integrated over must be piecewise linear: h(u2 | u1) is
  # the integral over t of the density at (u1, t), whatever udp1 is
  base <- bd("gumbel", 270, 2)
  N <- 20000
  for (x1 in list(v2p(0.4, 1.5), udplegendre(3), udpcosinebex(c(0.3, -0.8, 0.5)))) {
    o <- bsicopula(base, x1, udpzigzag(widths = c(3, 4, 3)))
    for (u1 in c(0.15, 0.6, 0.85)) {
      for (u2 in c(0.3, 0.55, 0.9)) {
        t <- (seq_len(N) - 0.5) / N * u2
        expect_equal(hbsicopula(u1, u2, o), mean(dbsicopula(rep(u1, N), t, o)) * u2, tolerance = 1e-3)
      }
    }
    # conversely, conditioning on U2 integrates over the non-linear udp1
    t <- (seq_len(N) - 0.5) / N * 0.55
    expect_equal(hbsicopula(0.55, 0.4, o, cond_var = 2), mean(dbsicopula(t, rep(0.4, N), o)) * 0.55, tolerance = 1e-3)
  }
})

test_that("a linear v-transform commutes with the h-function: T(h) = h*(T, T)", {
  skip_if_not_installed("rvinecopulib")
  set.seed(6)
  u <- runif(300)
  v <- runif(300)
  b <- bd("gumbel", 270, 2)
  for (T in list(vlinear(0.4), udpzigzag(widths = c(3, 4, 3)), udpzigzag(widths = c(1, 2, 1), up = TRUE), udpcosine(4))) {
    o <- bsicopula(b, T, T)
    Tu <- udptrans(T, u)
    Tv <- udptrans(T, v)
    expect_equal(udptrans(T, hbsicopula(u, v, o, cond_var = 1)), rvinecopulib::hbicop(cbind(Tu, Tv), 1, b), tolerance = 1e-8)
    expect_equal(udptrans(T, hbsicopula(u, v, o, cond_var = 2)), rvinecopulib::hbicop(cbind(Tu, Tv), 2, b), tolerance = 1e-8)
  }
})

test_that("a parCopula base copula agrees with bicop_dist, for both conditioning variables", {
  skip_if_not_installed("rvinecopulib")
  skip_if_not_installed("copula")
  set.seed(7)
  u <- runif(60)
  v <- runif(60)
  pairs <- list(
    list(copula::claytonCopula(2), bd("clayton", 0, 2)),
    list(copula::normalCopula(0.5), bd("gaussian", 0, 0.5)),
    list(copula::gumbelCopula(3), bd("gumbel", 0, 3))
  )
  for (pr in pairs) {
    a <- bsicopula(pr[[1]], vlinear(0.4), udpzigzag(widths = c(3, 4, 3)))
    b <- bsicopula(pr[[2]], vlinear(0.4), udpzigzag(widths = c(3, 4, 3)))
    expect_equal(pbsicopula(u, v, a), pbsicopula(u, v, b), tolerance = 1e-7)
    expect_equal(hbsicopula(u, v, a), hbsicopula(u, v, b), tolerance = 1e-6)
    expect_equal(hbsicopula(u, v, a, cond_var = 2), hbsicopula(u, v, b, cond_var = 2), tolerance = 1e-6)
    expect_equal(hbsicopula(u, v, a, inverse = TRUE), hbsicopula(u, v, b, inverse = TRUE), tolerance = 1e-6)
    expect_equal(hbsicopula(u, v, a, cond_var = 2, inverse = TRUE),
      hbsicopula(u, v, b, cond_var = 2, inverse = TRUE),
      tolerance = 1e-6
    )
  }
})

test_that("a non-exchangeable parCopula is handled consistently in both directions", {
  skip_if_not_installed("rvinecopulib")
  skip_if_not_installed("copula")
  # a rotation that flips the second coordinate: copula::cCopula() would give
  # the complement of the conditional distribution here, so the h-functions
  # are checked against the derivative of the CDF and against simulation
  set.seed(8)
  for (flip in list(c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))) {
    o <- bsicopula(copula::rotCopula(copula::claytonCopula(2), flip = flip), vlinear(0.4), udpcosine(3))
    u <- runif(80, 0.05, 0.95)
    v <- runif(80, 0.05, 0.95)
    d <- 1e-5
    expect_equal(hbsicopula(u, v, o), (pbsicopula(u + d, v, o) - pbsicopula(u - d, v, o)) / (2 * d), tolerance = 1e-4)
    expect_equal(hbsicopula(u, v, o, cond_var = 2), (pbsicopula(u, v + d, o) - pbsicopula(u, v - d, o)) / (2 * d), tolerance = 1e-4)
    p <- runif(80)
    expect_equal(hbsicopula(u, hbsicopula(u, p, o, inverse = TRUE), o), p, tolerance = 1e-6)
    expect_equal(hbsicopula(hbsicopula(p, v, o, cond_var = 2, inverse = TRUE), v, o, cond_var = 2), p, tolerance = 1e-6)
    U <- rbsicopula(2e5, o)
    sim <- mean(U[, 1] <= 0.4 & U[, 2] <= 0.6)
    expect_lt(abs(pbsicopula(0.4, 0.6, o) - sim), 0.004)
  }
})

test_that("arguments are validated and unsupported models rejected", {
  skip_if_not_installed("rvinecopulib")
  o <- linear_models()[[1]]
  expect_error(pbsicopula(0.3, 0.4, "x"), "class 'bsicopula'")
  expect_error(pbsicopula(0.3, 0.4), "class 'bsicopula'")
  expect_error(hbsicopula(0.3, 0.4, o, cond_var = 3), "'cond_var' must be 1 or 2")
  expect_error(hbsicopula(0.3, 0.4, o, inverse = NA), "'inverse' must be TRUE or FALSE")
  expect_error(pbsicopula(1.2, 0.4, o), "numbers in \\[0, 1\\]")
  expect_error(pbsicopula(c(0.2, NA), c(0.4, 0.5), o), "numbers in \\[0, 1\\]")
  expect_error(pbsicopula(c(0.1, 0.2, 0.3), c(0.4, 0.5), o), "same length")
  expect_error(pbsicopula(0.3, 0.4, o, nodes = 2), "'nodes' must be")
  expect_error(pbsicopula(c(0.1, 0.2, 0.3), object = o), "two-column matrix")
  # non-linear margins are handled by quadrature
  nl <- bsicopula(bd("gaussian", 0, 0.5), vlinear(0.4), v2p(0.4, 1.5))
  expect_silent(pbsicopula(0.3, 0.4, nl))
  expect_silent(hbsicopula(0.3, 0.4, nl))
})

test_that("input may be a matrix, with the object in the second position, and length-1 arguments recycle", {
  skip_if_not_installed("rvinecopulib")
  o <- linear_models()[[1]]
  set.seed(9)
  U <- cbind(runif(10), runif(10))
  expect_equal(pbsicopula(U, object = o), pbsicopula(U[, 1], U[, 2], o))
  expect_equal(pbsicopula(U, o), pbsicopula(U[, 1], U[, 2], o))
  expect_equal(pbsicopula(as.data.frame(U), o), pbsicopula(U[, 1], U[, 2], o))
  expect_equal(hbsicopula(U, object = o, cond_var = 2, inverse = TRUE), hbsicopula(U[, 1], U[, 2], o, cond_var = 2, inverse = TRUE))
  expect_equal(dbsicopula(U, o), dbsicopula(U[, 1], U[, 2], o))
  expect_equal(dbsicopula(U, object = o), dbsicopula(U[, 1], U[, 2], o))
  g <- seq(0.1, 0.9, 0.2)
  expect_equal(hbsicopula(0.4, g, o), hbsicopula(rep(0.4, 5), g, o))
  expect_equal(pbsicopula(g, 0.4, o), pbsicopula(g, rep(0.4, 5), o))
  expect_equal(dbsicopula(g, 0.4, o), dbsicopula(g, rep(0.4, 5), o))
  expect_identical(length(pbsicopula(numeric(0), numeric(0), o)), 0L)
  expect_identical(length(hbsicopula(numeric(0), numeric(0), o, inverse = TRUE)), 0L)
  expect_error(dbsicopula(c(0.1, 0.2, 0.3), c(0.4, 0.5), o), "same length")
})

# ---- non-linear transformations (tanh-sinh quadrature) ---------------------

nonlinear_models <- function() {
  list(
    `gumbel270 | v2p x v3p` = bsicopula(bd("gumbel", 270, 2), v2p(0.4, 1.5), v3p(0.6, 1.3, 1.2)),
    `clayton90 | vlinear x v2b` = bsicopula(bd("clayton", 90, 2), vlinear(0.4), v2b(0.3, 1.4)),
    `gaussian | legendre x v2p` = bsicopula(bd("gaussian", 0, 0.6), udplegendre(3), v2p(0.5, 0.8))
  )
}

test_that("quadrature h agrees with integrate() of the density", {
  skip_if_not_installed("rvinecopulib")
  set.seed(2)
  u1 <- runif(30)
  u2 <- runif(30)
  for (nm in names(nonlinear_models())) {
    bc <- nonlinear_models()[[nm]]
    h <- hbsicopula(u1, u2, bc)
    br <- udp:::udpbreaks(bc@udp2)
    ref <- vapply(seq_along(u1), function(i) {
      v1 <- udptrans(bc@udp1, u1[i])
      g <- function(x) rvinecopulib::dbicop(cbind(rep(v1, length(x)), udptrans(bc@udp2, x)), bc@basecopula)
      sum(vapply(seq_len(length(br) - 1L), function(k) {
        hi <- min(max(u2[i], br[k]), br[k + 1L])
        if (hi > br[k]) integrate(g, br[k], hi, rel.tol = 1e-10, subdivisions = 1000L)$value else 0
      }, 0))
    }, 0)
    expect_lt(max(abs(h - ref)), 1e-5, label = nm)
  }
})

test_that("quadrature and exact tiers agree on a piecewise-linear target", {
  skip_if_not_installed("rvinecopulib")
  set.seed(3)
  bc <- bsicopula(bd("clayton", 90, 2), vlinear(0.4), udpzigzag(widths = c(3, 4, 3)))
  u1 <- runif(100)
  u2 <- runif(100)
  q <- udp:::h_quad(u2, udptrans(bc@udp1, u1), bc@udp2, bc, 1L, 101L)
  expect_lt(max(abs(q - hbsicopula(u1, u2, bc))), 1e-7)
})

test_that("inverse h inverts h for non-linear transformations, both conditioning variables", {
  skip_if_not_installed("rvinecopulib")
  set.seed(4)
  u1 <- runif(50)
  u2 <- runif(50)
  for (bc in nonlinear_models()) {
    for (cv in 1:2) {
      h <- hbsicopula(u1, u2, bc, cond_var = cv)
      expect_lt(max(abs(hbsicopula(u1, u2, bc, cond_var = cv) - h)), 1e-15)
      back <- if (cv == 1L) hbsicopula(u1, h, bc, cond_var = 1, inverse = TRUE) else hbsicopula(h, u2, bc, cond_var = 2, inverse = TRUE)
      target <- if (cv == 1L) u2 else u1
      expect_lt(max(abs(back - target)), 1e-6)
    }
  }
})

test_that("non-linear CDF matches simulation, boundaries and the derivative relation", {
  skip_if_not_installed("rvinecopulib")
  set.seed(5)
  pts <- cbind(c(0.2, 0.5, 0.8, 0.4), c(0.3, 0.5, 0.6, 0.9))
  for (bc in nonlinear_models()) {
    sim <- rbsicopula(2e5, bc)
    est <- vapply(1:4, function(i) mean(sim[, 1] <= pts[i, 1] & sim[, 2] <= pts[i, 2]), 0)
    expect_lt(max(abs(pbsicopula(pts[, 1], pts[, 2], bc) - est)), 4e-3)
    expect_equal(pbsicopula(c(0, 1, 1, 0.5), c(0.5, 0.5, 1, 1), bc), c(0, 0.5, 1, 0.5), tolerance = 1e-6)
    e <- 1e-4
    d <- (pbsicopula(0.5 + e, 0.6, bc) - pbsicopula(0.5 - e, 0.6, bc)) / (2 * e)
    expect_equal(d, hbsicopula(0.5, 0.6, bc), tolerance = 1e-5)
  }
})

# ---- randomizer models ------------------------------------------------------

randomizer_models <- function() {
  rs <- randsdvine(bd("gaussian", 0, 0.7), bd("clayton", 90, 1.5), bd("gumbel", 270, 1.8))
  list(
    `randsdvine | v2p x v3p` = bsicopula(bd("gumbel", 270, 2), v2p(0.4, 1.5), v3p(0.6, 1.3, 1.2), rs),
    `randsdvine | vlinear x zigzag` = bsicopula(bd("clayton", 90, 2), vlinear(0.4), udpzigzag(widths = c(3, 4, 3)), rs),
    `randmixture | vlinear x cosine` = bsicopula(
      bd("gaussian", 0, 0.5), vlinear(0.3), udpcosine(3),
      randmixture(bd("gaussian", 0, 0.6), bd("clayton", 0, 2), function(v1, v2) v1 + v2 < 1)
    )
  )
}

test_that("h with a randomizer integrates dbsicopula(), whichever margin is conditioned on", {
  skip_if_not_installed("rvinecopulib")
  set.seed(7)
  N <- 20000
  g <- (seq_len(N) - 0.5) / N
  for (nm in names(randomizer_models())) {
    bc <- randomizer_models()[[nm]]
    tol <- if (grepl("mixture", nm)) 2e-3 else 1e-4 # the mixture selector is discontinuous in v
    for (i in 1:4) {
      u <- runif(1)
      x <- runif(1)
      expect_equal(hbsicopula(u, x, bc), mean(dbsicopula(rep(u, N), g * x, bc)) * x, tolerance = tol, label = nm)
      expect_equal(hbsicopula(x, u, bc, cond_var = 2), mean(dbsicopula(g * x, rep(u, N), bc)) * x, tolerance = tol, label = nm)
    }
    # uniform margins: h(1 | u) = 1 and C(u, 1) = u
    u <- runif(5)
    expect_equal(hbsicopula(u, rep(1, 5), bc), rep(1, 5), tolerance = 1e-6)
    expect_equal(pbsicopula(u, rep(1, 5), bc, nodes = 41), u, tolerance = 2e-3)
  }
})

test_that("h with a randomizer has an inverse (in the sense of h), and CDF matches simulation", {
  skip_if_not_installed("rvinecopulib")
  set.seed(8)
  u1 <- runif(30)
  u2 <- runif(30)
  pts <- cbind(c(0.2, 0.5, 0.8, 0.4), c(0.3, 0.5, 0.6, 0.9))
  for (nm in names(randomizer_models())) {
    bc <- randomizer_models()[[nm]]
    tol <- if (grepl("mixture", nm)) 3e-3 else 1e-4
    # h can be flat where the density vanishes, so compare h values, not arguments
    h <- hbsicopula(u1, u2, bc)
    expect_lt(max(abs(hbsicopula(u1, hbsicopula(u1, h, bc, inverse = TRUE), bc) - h)), tol, label = nm)
    h2 <- hbsicopula(u1, u2, bc, cond_var = 2)
    back <- hbsicopula(h2, u2, bc, cond_var = 2, inverse = TRUE)
    expect_lt(max(abs(hbsicopula(back, u2, bc, cond_var = 2) - h2)), tol, label = nm)
    sim <- rbsicopula(2e5, bc)
    est <- vapply(1:4, function(i) mean(sim[, 1] <= pts[i, 1] & sim[, 2] <= pts[i, 2]), 0)
    expect_lt(max(abs(pbsicopula(pts[, 1], pts[, 2], bc, nodes = 41) - est)), 4e-3, label = nm)
  }
})
