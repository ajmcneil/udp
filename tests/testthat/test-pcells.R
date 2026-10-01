# pcells(): cell probabilities of a bsicopula. Needs rvinecopulib throughout.

pcells_models <- function() {
  bd <- rvinecopulib::bicop_dist
  rm3 <- randsdvine(bd("gaussian", 0, 0.7), bd("clayton", 0, 1.5), bd("gumbel", 0, 1.8))
  list(
    vv = bsicopula(bd("gaussian", 0, 0.5), vlinear(0.5), vlinear(0.5), rm3),
    zz = bsicopula(bd("clayton", 90, 2), vlinear(0.4), udpzigzag(widths = c(0.3, 0.45, 0.25)), rm3),
    lb = bsicopula(bd("gaussian", 0, 0.5), udpcosine(2), udplegendrebex(c(0.3, -0.8, 0.5)), rm3),
    l4 = bsicopula(bd("gumbel", 0, 2), udplegendre(4), udpcosinebex(c(0.3, -0.8, 0.5)), rm3),
    null = bsicopula(bd("gaussian", 0, 0.5), v2p(0.4, 1.5), udplegendrebex(c(0.3, -0.8, 0.5))),
    mix = bsicopula(bd("clayton", 0, 2), v2p(0.4, 1.5), udpzigzag(widths = c(0.3, 0.45, 0.25)),
      randmixture(bd("gaussian", 0, 0.8), bd("gaussian", 0, -0.8), function(v1, v2) pmax(v1, v2) > 0.6)
    )
  )
}

test_that("pcells() validates its arguments", {
  skip_if_not_installed("rvinecopulib")
  o <- pcells_models()$vv
  expect_error(pcells("x"), "class 'bsicopula'")
  expect_error(pcells(o, 1.2, 0.5), "'v1' must be NULL or a numeric vector")
  expect_error(pcells(o, 0.2, NA_real_), "'v2' must be NULL or a numeric vector")
  expect_error(pcells(o, 0.2, 0.5, cells = "other"), "should be one of")
  expect_error(pcells(o, 0.2, 0.5, grid = NA), "'grid' must be TRUE or FALSE")
  expect_error(pcells(o, ngrid = 1), "'ngrid' must be")
})

test_that("conditional cell probabilities sum to 1, including at turning values", {
  skip_if_not_installed("rvinecopulib")
  set.seed(1)
  v1 <- c(0, 1, 0, 1, runif(20))
  v2 <- c(0, 0, 1, 1, runif(20))
  for (o in pcells_models()) {
    for (cells in c("monotone", "smooth")) {
      P <- pcells(o, v1, v2, cells = cells)
      expect_equal(dim(P)[3], length(v1))
      expect_true(all(P >= -1e-12))
      expect_equal(apply(P, 3, sum), rep(1, length(v1)), tolerance = 1e-12)
    }
  }
})

test_that("with independent randomizers the probabilities are products of the marginal ones", {
  skip_if_not_installed("rvinecopulib")
  o <- pcells_models()$null
  b1 <- udpmonobreaks(o@udp1)
  b2 <- udpmonobreaks(o@udp2)
  marginal <- function(x, v, breaks) {
    M <- udpinverse(x, v, prob = TRUE)
    keep <- !is.na(M[1, ])
    cell <- findInterval(M[1, keep], breaks, all.inside = TRUE)
    vapply(seq_len(length(breaks) - 1L), function(l) sum(attr(M, "prob")[1, keep][cell == l]), 0)
  }
  for (v in list(c(0.3, 0.6), c(0.85, 0.12))) {
    expect_equal(
      unname(unclass(pcells(o, v[1], v[2]))[, , drop = TRUE]),
      outer(marginal(o@udp1, v[1], b1), marginal(o@udp2, v[2], b2)),
      ignore_attr = TRUE
    )
  }
})

test_that("at a v-transform's fulcrum value the branches keep their own cells", {
  skip_if_not_installed("rvinecopulib")
  o <- pcells_models()$vv
  # both pre-images sit on the fulcrum at v = 0: the down branch is cell 1
  Fz <- randomizer_cdf(0, 0, o@basecopula, o@randomizermod)
  expect_equal(pcells(o, 0, 0)[1, 1], Fz(0.5, 0.5))
  expect_equal(sum(pcells(o, 0, 0)), 1)
})

test_that("cells = 'smooth' refines cells = 'monotone'", {
  skip_if_not_installed("rvinecopulib")
  o <- pcells_models()$l4
  Pm <- pcells(o, 0.3, 0.7)
  Ps <- pcells(o, 0.3, 0.7, cells = "smooth")
  expect_identical(dim(Pm), c(4L, 3L))
  expect_true(all(dim(Ps) > dim(Pm)))
  # adding up the smooth cells inside each monotone cell recovers the latter
  within <- function(fine, coarse) {
    findInterval((fine[-1] + fine[-length(fine)]) / 2, coarse, all.inside = TRUE)
  }
  g1 <- within(attr(Ps, "breaks1"), attr(Pm, "breaks1"))
  g2 <- within(attr(Ps, "breaks2"), attr(Pm, "breaks2"))
  agg <- rowsum(t(rowsum(unclass(Ps)[, ], g1)), g2)
  expect_equal(t(agg), unclass(Pm)[, ], ignore_attr = TRUE)
  # the two coincide for classes without extra break points
  z <- pcells_models()$zz
  expect_equal(pcells(z, 0.3, 0.7), pcells(z, 0.3, 0.7, cells = "smooth"))
})

test_that("grid = TRUE evaluates every combination, and drop controls the shape", {
  skip_if_not_installed("rvinecopulib")
  o <- pcells_models()$zz
  v1 <- c(0.1, 0.4, 0.8)
  v2 <- c(0.25, 0.6)
  G <- pcells(o, v1, v2, grid = TRUE)
  expect_identical(dim(G), c(2L, 3L, 3L, 2L))
  for (i in seq_along(v1)) {
    for (j in seq_along(v2)) {
      expect_equal(G[, , i, j], pcells(o, v1[i], v2[j])[, ], ignore_attr = TRUE)
    }
  }
  expect_identical(dim(pcells(o, 0.3, 0.6)), c(2L, 3L))
  expect_identical(dim(pcells(o, 0.3, 0.6, drop = FALSE)), c(2L, 3L, 1L))
  expect_identical(dim(pcells(o, v1, 0.6)), c(2L, 3L, 3L))
  expect_identical(dim(pcells(o, 0.3, 0.6, grid = TRUE, drop = FALSE)), c(2L, 3L, 1L, 1L))
  expect_identical(dim(pcells(o, ngrid = 20)), c(2L, 3L))
  expect_identical(dim(pcells(o, ngrid = 20, drop = FALSE)), c(2L, 3L, 1L))
  expect_identical(dim(pcells(o, v2 = v2, ngrid = 20)), c(2L, 3L, 2L))
  expect_identical(attr(G, "breaks2"), c(0, 0.3, 0.75, 1))
})

test_that("unconditional probabilities have the cell widths as row and column sums", {
  skip_if_not_installed("rvinecopulib")
  for (o in pcells_models()) {
    P <- pcells(o, ngrid = 200)
    expect_equal(sum(P), 1, tolerance = 1e-10)
    expect_equal(unname(rowSums(P)), diff(attr(P, "breaks1")), tolerance = 2e-3)
    expect_equal(unname(colSums(P)), diff(attr(P, "breaks2")), tolerance = 2e-3)
  }
})

test_that("conditioning on one carrier averages the fully conditional probabilities", {
  skip_if_not_installed("rvinecopulib")
  o <- pcells_models()$zz
  # given V1: the average over V2 | V1 = v1, here by brute-force conditional quantiles
  w <- (seq_len(2000) - 0.5) / 2000
  V2 <- rvinecopulib::hbicop(cbind(0.3, w), cond_var = 1, family = o@basecopula, inverse = TRUE)
  brute <- rowMeans(pcells(o, 0.3, V2), dims = 2L)
  expect_equal(unclass(pcells(o, v1 = 0.3))[, ], brute, tolerance = 1e-3, ignore_attr = TRUE)
  V1 <- rvinecopulib::hbicop(cbind(w, 0.7), cond_var = 2, family = o@basecopula, inverse = TRUE)
  brute <- rowMeans(pcells(o, V1, 0.7), dims = 2L)
  expect_equal(unclass(pcells(o, v2 = 0.7))[, ], brute, tolerance = 1e-3, ignore_attr = TRUE)
  # each slice still sums to 1, and the margin that was conditioned on keeps
  # its own selection probabilities
  P <- pcells(o, v1 = c(0.2, 0.6))
  expect_equal(apply(P, 3, sum), c(1, 1), tolerance = 1e-10)
  expect_equal(unname(rowSums(P[, , 1])), c(0.4, 0.6), tolerance = 1e-6)
})

test_that("conditional probabilities agree with simulated selections, for unsymmetric copulas in every slot", {
  skip_if_not_installed("rvinecopulib")
  # regression: the randsdvine code once evaluated P(V1 <= v1 | V2 = v2) with
  # the base copula's arguments swapped, and used copZ1V2_V1 in the order
  # (V2, Z1) instead of (Z1, V2) -- both invisible for exchangeable copulas
  bd <- rvinecopulib::bicop_dist
  set.seed(2)
  N <- 1e5
  for (rm in list(
    randsdvine(bd("gaussian", 0, 0.7), bd("clayton", 0, 1.5), bd("gumbel", 0, 1.8)),
    randsdvine(bd("clayton", 90, 3), bd("clayton", 90, 3), bd("gumbel", 270, 2)),
    randsdvine(bd("clayton", 270, 3), bd("clayton", 270, 3), bd("clayton", 90, 2))
  )) {
    o <- bsicopula(bd("clayton", 90, 2), vlinear(0.4), udpzigzag(widths = c(0.3, 0.45, 0.25)), rm)
    for (v in list(c(0.3, 0.5), c(0.7, 0.2))) {
      Z <- randsdvine_sample(rep(v[1], N), rep(v[2], N), o@basecopula, o@randomizermod)
      U1 <- udpsi(o@udp1, rep(v[1], N), Z[, "Z1"])
      U2 <- udpsi(o@udp2, rep(v[2], N), Z[, "Z2"])
      P <- pcells(o, v[1], v[2])
      sim <- table(
        factor(findInterval(U1, attr(P, "breaks1"), all.inside = TRUE), 1:2),
        factor(findInterval(U2, attr(P, "breaks2"), all.inside = TRUE), 1:3)
      ) / N
      expect_lt(max(abs(unclass(P)[, ] - unclass(sim)[, ])), 0.006)
    }
  }
})

test_that("a parCopula base copula gives the same integrated probabilities as a bicop_dist", {
  skip_if_not_installed("rvinecopulib")
  skip_if_not_installed("copula")
  mk <- function(cop) bsicopula(cop, v2p(0.4, 1.5), udplegendrebex(c(0.3, -0.8, 0.5)))
  a <- mk(copula::claytonCopula(2))
  b <- mk(rvinecopulib::bicop_dist("clayton", 0, 2))
  expect_equal(pcells(a, 0.3, 0.6), pcells(b, 0.3, 0.6))
  expect_equal(pcells(a, v1 = c(0.1, 0.5)), pcells(b, v1 = c(0.1, 0.5)), tolerance = 5e-3)
  expect_equal(pcells(a, v2 = c(0.1, 0.5)), pcells(b, v2 = c(0.1, 0.5)), tolerance = 5e-3)
  expect_equal(pcells(a, ngrid = 200), pcells(b, ngrid = 200), tolerance = 5e-3)
})
