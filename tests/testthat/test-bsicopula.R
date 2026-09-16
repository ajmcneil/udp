# sdvine(), bsicopula(), rbsicopula() -- these need the Suggested packages
# copula and rvinecopulib, so every test here skips cleanly when either is
# unavailable.

test_that("sdvine() requires bicop_dist objects and defaults the two tree-2 edges to independence", {
  skip_if_not_installed("rvinecopulib")
  cop3 <- rvinecopulib::bicop_dist("t", 0, c(0.4, 5))
  sdv <- sdvine(cop3)
  expect_s4_class(sdv, "sdvine")
  expect_identical(sdv@copZ1Z2_V1V2$family, "t")
  expect_identical(sdv@copZ1V2_V1$family, "indep")
  expect_identical(sdv@copV1Z2_V2$family, "indep")

  expect_error(sdvine(rvinecopulib::bicop_dist("clayton", 0, 2), copZ1V2_V1 = 1),
    "bicop_dist"
  )
  skip_if_not_installed("copula")
  expect_error(sdvine(copula::claytonCopula(2)), "bicop_dist")
})

test_that("bsicopula() accepts either a parCopula or bicop_dist basecopula when randomizermod is NULL", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  bc1 <- bsicopula(copula::claytonCopula(2), udpcosine(2), udpcosine(3))
  bc2 <- bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3))
  expect_s4_class(bc1, "bsicopula")
  expect_s4_class(bc2, "bsicopula")
  expect_null(bc1@randomizermod)
})

test_that("bsicopula() rejects non-udp margins and non-copula basecopula", {
  skip_if_not_installed("rvinecopulib")
  bic <- rvinecopulib::bicop_dist("clayton", 0, 2)
  expect_error(bsicopula(bic, 1, udpcosine(3)), "class 'udp'")
  expect_error(bsicopula(bic, udpcosine(2), "not a udp"), "class 'udp'")
  expect_error(bsicopula(matrix(1), udpcosine(2), udpcosine(3)), "parCopula.*bicop_dist")
})

test_that("bsicopula() requires basecopula to be bicop_dist when randomizermod is an sdvine", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  sdv <- sdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)))
  expect_error(
    bsicopula(copula::claytonCopula(2), udpcosine(2), udpcosine(3), randomizermod = sdv),
    "bicop_dist object when 'randomizermod'"
  )
  expect_error(
    bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3), randomizermod = "x"),
    "class 'sdvine'"
  )
  expect_s4_class(
    bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3), randomizermod = sdv),
    "bsicopula"
  )
})

test_that("rbsicopula() requires a bsicopula object", {
  expect_error(rbsicopula(10, list()), "class 'bsicopula'")
})

test_that("rbsicopula() with randomizermod = NULL returns an n x 2 matrix, uniform margins", {
  skip_if_not_installed("rvinecopulib")
  set.seed(1)
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3))
  samp <- rbsicopula(5000, bc)
  expect_equal(dim(samp), c(5000L, 2L))
  expect_identical(colnames(samp), c("U1", "U2"))
  expect_true(all(samp >= 0 & samp <= 1))
  expect_gt(ks.test(samp[, "U1"], "punif")$p.value, 0.001)
  expect_gt(ks.test(samp[, "U2"], "punif")$p.value, 0.001)
})

test_that("rbsicopula() with an sdvine matches a direct hand-rolled implementation of the algorithm", {
  skip_if_not_installed("rvinecopulib")
  bic <- rvinecopulib::bicop_dist
  hbicop <- rvinecopulib::hbicop
  basecop <- bic("gumbel", 0, 1.6)
  sdv <- sdvine(bic("t", 0, c(0.4, 5)), bic("clayton", 0, 1.3), bic("joe", 0, 2.0))

  set.seed(11)
  n <- 500
  V1 <- runif(n)
  V2 <- runif(n)

  set.seed(22)
  w1 <- runif(n)
  w2 <- runif(n)
  e21 <- hbicop(cbind(V1, V2), cond_var = 1, family = basecop)
  e12 <- hbicop(cbind(V2, V1), cond_var = 1, family = basecop)
  z1 <- hbicop(cbind(e21, w1), cond_var = 1, family = sdv@copZ1V2_V1, inverse = TRUE)
  y <- hbicop(cbind(w1, w2), cond_var = 1, family = sdv@copZ1Z2_V1V2, inverse = TRUE)
  z2 <- hbicop(cbind(e12, y), cond_var = 1, family = sdv@copV1Z2_V2, inverse = TRUE)

  set.seed(22)
  Z <- sdvine_sample(V1, V2, basecop, sdv)

  expect_equal(as.numeric(Z[, "Z1"]), z1)
  expect_equal(as.numeric(Z[, "Z2"]), z2)
})

test_that("rbsicopula() with an sdvine still gives uniform margins and matches udpsi() directly", {
  skip_if_not_installed("rvinecopulib")
  sdv <- sdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)),
    rvinecopulib::bicop_dist("clayton", 0, 1.3),
    rvinecopulib::bicop_dist("joe", 0, 2.0)
  )
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3),
    randomizermod = sdv
  )

  set.seed(3)
  samp <- rbsicopula(5000, bc)
  expect_equal(dim(samp), c(5000L, 2L))
  expect_true(all(samp >= 0 & samp <= 1))
  expect_gt(ks.test(samp[, "U1"], "punif")$p.value, 0.001)
  expect_gt(ks.test(samp[, "U2"], "punif")$p.value, 0.001)

  set.seed(3)
  V <- rvinecopulib::rbicop(5000, bc@basecopula)
  Z <- sdvine_sample(V[, 1], V[, 2], bc@basecopula, sdv)
  expect_equal(samp[, "U1"], udpsi(udpcosine(2), V[, 1], Z[, "Z1"]))
  expect_equal(samp[, "U2"], udpsi(udpcosine(3), V[, 2], Z[, "Z2"]))
})
