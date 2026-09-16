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

test_that("dbsicopula() requires a bsicopula object and matching-length u1/u2", {
  skip_if_not_installed("rvinecopulib")
  bc <- bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3))
  expect_error(dbsicopula(0.2, 0.3, list()), "class 'bsicopula'")
  expect_error(dbsicopula(c(0.2, 0.3), 0.3, bc), "same length")
})

test_that("dbsicopula() with randomizermod = NULL is exactly c_V(T1(u1), T2(u2))", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  u1 <- c(0.1, 0.37, 0.6, 0.85)
  u2 <- c(0.2, 0.55, 0.4, 0.9)

  bc_par <- bsicopula(copula::claytonCopula(2), udpcosine(2), udpcosine(3))
  v1 <- udptrans(udpcosine(2), u1)
  v2 <- udptrans(udpcosine(3), u2)
  expect_equal(dbsicopula(u1, u2, bc_par), copula::dCopula(cbind(v1, v2), copula::claytonCopula(2)))

  bc_bicop <- bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3))
  expect_equal(
    dbsicopula(u1, u2, bc_bicop),
    rvinecopulib::dbicop(cbind(v1, v2), rvinecopulib::bicop_dist("clayton", 0, 2))
  )
})

test_that("dbsicopula() with an all-independence sdvine matches randomizermod = NULL", {
  skip_if_not_installed("rvinecopulib")
  u1 <- c(0.1, 0.37, 0.6, 0.85)
  u2 <- c(0.2, 0.55, 0.4, 0.9)
  base <- rvinecopulib::bicop_dist("gumbel", 0, 1.6)
  sdv_indep <- sdvine(rvinecopulib::bicop_dist(), rvinecopulib::bicop_dist(), rvinecopulib::bicop_dist())

  bc_indep <- bsicopula(base, udpcosine(2), udpcosine(3))
  bc_sdv <- bsicopula(base, udpcosine(2), udpcosine(3), randomizermod = sdv_indep)
  expect_equal(dbsicopula(u1, u2, bc_sdv), dbsicopula(u1, u2, bc_indep), tolerance = 1e-6)
})

test_that("dbsicopula() with a non-trivial sdvine integrates to 1 over the unit square", {
  skip_if_not_installed("rvinecopulib")
  sdv <- sdvine(
    rvinecopulib::bicop_dist("t", 0, c(0.4, 5)),
    rvinecopulib::bicop_dist("clayton", 0, 1.3),
    rvinecopulib::bicop_dist("joe", 0, 2.0)
  )
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3),
    randomizermod = sdv
  )
  # midpoint rule on the product grid of each margin's breakpoints, fine
  # enough to keep quadrature error well under the sampling noise this test
  # isn't otherwise trying to control for
  bp1 <- c(0, 1 / 2, 1)
  bp2 <- c(0, 1 / 3, 2 / 3, 1)
  ng <- 60
  total <- 0
  for (i in seq_len(length(bp1) - 1)) {
    for (j in seq_len(length(bp2) - 1)) {
      x1 <- bp1[i] + (bp1[i + 1] - bp1[i]) * ((seq_len(ng) - 0.5) / ng)
      x2 <- bp2[j] + (bp2[j + 1] - bp2[j]) * ((seq_len(ng) - 0.5) / ng)
      grid <- expand.grid(x1 = x1, x2 = x2)
      total <- total + mean(dbsicopula(grid$x1, grid$x2, bc)) * (bp1[i + 1] - bp1[i]) * (bp2[j + 1] - bp2[j])
    }
  }
  expect_equal(total, 1, tolerance = 0.01)
})

test_that("dbsicopula() gives a finite, non-negative value when u1 lands exactly on a breakpoint", {
  skip_if_not_installed("rvinecopulib")
  sdv <- sdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)))
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3),
    randomizermod = sdv
  )
  val <- dbsicopula(1 / 2, 0.4, bc)
  expect_true(is.finite(val) && val >= 0)
})

## plot() ---------------------------------------------------------------

test_that("plot() runs for both type = \"contour\" (default) and \"persp\", with and without an sdvine", {
  skip_if_not_installed("rvinecopulib")
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })

  sdv <- sdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)), rvinecopulib::bicop_dist("clayton", 0, 1.3))
  bc_indep <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3))
  bc_sdv <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3),
    randomizermod = sdv
  )
  for (bc in list(bc_indep, bc_sdv)) {
    expect_no_error(plot(bc))
    expect_no_error(plot(bc, type = "persp"))
    expect_no_error(plot(bc, n = 30))
  }
})

test_that("plot() rejects an unknown type", {
  skip_if_not_installed("rvinecopulib")
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3))
  expect_error(plot(bc, type = "nonsense"), "should be one of")
})

test_that("plot(type = \"contour\") defaults to quantile-based levels and honours a user-supplied 'levels'", {
  skip_if_not_installed("rvinecopulib")
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3))

  captured <- NULL
  local_mocked_bindings(
    contour = function(x, y, z, ..., levels) {
      captured <<- levels
      invisible(NULL)
    },
    .package = "udp"
  )
  plot(bc, n = 20)
  expect_length(captured, 9L) # deciles 0.1, ..., 0.9
  expect_true(!is.unsorted(captured, strictly = TRUE))

  plot(bc, n = 20, levels = c(0.5, 1, 2))
  expect_identical(captured, c(0.5, 1, 2))
})
