# randsdvine(), bsicopula(), rbsicopula() -- these need the Suggested packages
# copula and rvinecopulib, so every test here skips cleanly when either is
# unavailable.

test_that("randsdvine() requires bicop_dist objects and defaults the two tree-2 edges to independence", {
  skip_if_not_installed("rvinecopulib")
  cop3 <- rvinecopulib::bicop_dist("t", 0, c(0.4, 5))
  sdv <- randsdvine(cop3)
  expect_s4_class(sdv, "randsdvine")
  expect_identical(sdv@copZ1Z2_V1V2$family, "t")
  expect_identical(sdv@copZ1V2_V1$family, "indep")
  expect_identical(sdv@copV1Z2_V2$family, "indep")

  expect_error(randsdvine(rvinecopulib::bicop_dist("clayton", 0, 2), copZ1V2_V1 = 1),
    "bicop_dist"
  )
  skip_if_not_installed("copula")
  expect_error(randsdvine(copula::claytonCopula(2)), "bicop_dist")
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

test_that("bsicopula() requires basecopula to be bicop_dist when randomizermod is a randsdvine", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  sdv <- randsdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)))
  expect_error(
    bsicopula(copula::claytonCopula(2), udpcosine(2), udpcosine(3), randomizermod = sdv),
    "bicop_dist object when 'randomizermod'"
  )
  expect_error(
    bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3), randomizermod = "x"),
    "class 'randsdvine'"
  )
  expect_s4_class(
    bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3), randomizermod = sdv),
    "bsicopula"
  )
})

## randmixture() / randmixture_sample() --------------------------------------------

test_that("randmixture() validates cop1, cop2 and selector", {
  skip_if_not_installed("rvinecopulib")
  cop1 <- rvinecopulib::bicop_dist("gaussian", parameters = 1)
  cop2 <- rvinecopulib::bicop_dist("gaussian", parameters = -1)
  sel <- function(v1, v2) pmax(v1, v2) > 0.7

  vm <- randmixture(cop1, cop2, sel)
  expect_s4_class(vm, "randmixture")

  expect_error(randmixture(cop1, "not a copula", sel), "parCopula.*bicop_dist")
  expect_error(randmixture("not a copula", cop2, sel), "parCopula.*bicop_dist")
  expect_error(randmixture(cop1, cop2, "not a function"), "must be a function")
})

test_that("bsicopula() accepts a randmixture with either basecopula type (no bicop_dist restriction)", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  cop1 <- rvinecopulib::bicop_dist("gaussian", parameters = 1)
  cop2 <- rvinecopulib::bicop_dist("gaussian", parameters = -1)
  vm <- randmixture(cop1, cop2, function(v1, v2) pmax(v1, v2) > 0.7)

  expect_s4_class(
    bsicopula(copula::claytonCopula(2), udpcosine(2), udpcosine(3), randomizermod = vm),
    "bsicopula"
  )
  expect_s4_class(
    bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3), randomizermod = vm),
    "bsicopula"
  )
})

test_that("randmixture_sample() gives Z1 independent of V1 and Z2 independent of V2, for any selector", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  cop1 <- rvinecopulib::bicop_dist("gaussian", parameters = 1)
  cop2 <- rvinecopulib::bicop_dist("gaussian", parameters = -1)
  vm <- randmixture(cop1, cop2, function(v1, v2) pmax(v1, v2) > 0.7)

  set.seed(1)
  n <- 20000
  V <- copula::rCopula(n, copula::claytonCopula(2))
  Z <- randmixture_sample(V[, 1], V[, 2], vm)

  expect_equal(cor(Z[, "Z1"], V[, 1]), 0, tolerance = 0.02)
  expect_equal(cor(Z[, "Z2"], V[, 2]), 0, tolerance = 0.02)
  expect_gt(ks.test(Z[, "Z1"], "punif")$p.value, 0.001)
  expect_gt(ks.test(Z[, "Z2"], "punif")$p.value, 0.001)
})

test_that("randmixture_sample() actually switches regimes according to 'selector'", {
  skip_if_not_installed("copula")
  skip_if_not_installed("rvinecopulib")
  cop1 <- rvinecopulib::bicop_dist("gaussian", parameters = 1) # comonotonic
  cop2 <- rvinecopulib::bicop_dist("gaussian", parameters = -1) # countermonotonic
  vm <- randmixture(cop1, cop2, function(v1, v2) pmax(v1, v2) > 0.7)

  set.seed(1)
  n <- 20000
  V <- copula::rCopula(n, copula::claytonCopula(2))
  sel <- pmax(V[, 1], V[, 2]) > 0.7
  Z <- randmixture_sample(V[, 1], V[, 2], vm)

  expect_equal(cor(Z[sel, "Z1"], Z[sel, "Z2"]), 1, tolerance = 1e-8)
  expect_equal(cor(Z[!sel, "Z1"], Z[!sel, "Z2"]), -1, tolerance = 1e-8)
})

test_that("randmixture_sample() errors clearly if 'selector' returns the wrong shape", {
  skip_if_not_installed("rvinecopulib")
  cop1 <- rvinecopulib::bicop_dist("gaussian", parameters = 0.5)
  vm_wronglen <- randmixture(cop1, cop1, function(v1, v2) TRUE)
  vm_notlogical <- randmixture(cop1, cop1, function(v1, v2) v1)

  expect_error(
    randmixture_sample(runif(5), runif(5), vm_wronglen),
    "logical vector the same length"
  )
  expect_error(
    randmixture_sample(runif(5), runif(5), vm_notlogical),
    "logical vector the same length"
  )
})

test_that("rbsicopula() with a randmixture gives an n x 2 matrix with uniform margins", {
  skip_if_not_installed("rvinecopulib")
  cop1 <- rvinecopulib::bicop_dist("gaussian", parameters = 1)
  cop2 <- rvinecopulib::bicop_dist("gaussian", parameters = -1)
  vm <- randmixture(cop1, cop2, function(v1, v2) pmax(v1, v2) > 0.7)
  bc <- bsicopula(rvinecopulib::bicop_dist("clayton", 0, 2), udpcosine(2), udpcosine(3),
    randomizermod = vm
  )

  set.seed(1)
  samp <- rbsicopula(5000, bc)
  expect_equal(dim(samp), c(5000L, 2L))
  expect_true(all(samp >= 0 & samp <= 1))
  expect_gt(ks.test(samp[, "U1"], "punif")$p.value, 0.001)
  expect_gt(ks.test(samp[, "U2"], "punif")$p.value, 0.001)
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

test_that("rbsicopula() with a randsdvine matches a direct hand-rolled implementation of the algorithm", {
  skip_if_not_installed("rvinecopulib")
  bic <- rvinecopulib::bicop_dist
  hbicop <- rvinecopulib::hbicop
  basecop <- bic("gumbel", 0, 1.6)
  sdv <- randsdvine(bic("t", 0, c(0.4, 5)), bic("clayton", 0, 1.3), bic("joe", 0, 2.0))

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
  Z <- randsdvine_sample(V1, V2, basecop, sdv)

  expect_equal(as.numeric(Z[, "Z1"]), z1)
  expect_equal(as.numeric(Z[, "Z2"]), z2)
})

test_that("rbsicopula() with a randsdvine still gives uniform margins and matches udpsi() directly", {
  skip_if_not_installed("rvinecopulib")
  sdv <- randsdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)),
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
  Z <- randsdvine_sample(V[, 1], V[, 2], bc@basecopula, sdv)
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

test_that("dbsicopula() with an all-independence randsdvine matches randomizermod = NULL", {
  skip_if_not_installed("rvinecopulib")
  u1 <- c(0.1, 0.37, 0.6, 0.85)
  u2 <- c(0.2, 0.55, 0.4, 0.9)
  base <- rvinecopulib::bicop_dist("gumbel", 0, 1.6)
  sdv_indep <- randsdvine(rvinecopulib::bicop_dist(), rvinecopulib::bicop_dist(), rvinecopulib::bicop_dist())

  bc_indep <- bsicopula(base, udpcosine(2), udpcosine(3))
  bc_sdv <- bsicopula(base, udpcosine(2), udpcosine(3), randomizermod = sdv_indep)
  expect_equal(dbsicopula(u1, u2, bc_sdv), dbsicopula(u1, u2, bc_indep), tolerance = 1e-6)
})

test_that("dbsicopula() with a non-trivial randsdvine integrates to 1 over the unit square", {
  skip_if_not_installed("rvinecopulib")
  sdv <- randsdvine(
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

test_that("dbsicopula() with an all-independence randmixture matches randomizermod = NULL", {
  skip_if_not_installed("rvinecopulib")
  u1 <- c(0.1, 0.37, 0.6, 0.85)
  u2 <- c(0.2, 0.55, 0.4, 0.9)
  base <- rvinecopulib::bicop_dist("gumbel", 0, 1.6)
  vm_indep <- randmixture(rvinecopulib::bicop_dist(), rvinecopulib::bicop_dist(),
    selector = function(v1, v2) v1 > 0.5 # any selector: both branches are independence
  )

  bc_indep <- bsicopula(base, udpcosine(2), udpcosine(3))
  bc_vmix <- bsicopula(base, udpcosine(2), udpcosine(3), randomizermod = vm_indep)
  expect_equal(dbsicopula(u1, u2, bc_vmix), dbsicopula(u1, u2, bc_indep), tolerance = 1e-6)
})

test_that("dbsicopula() with a non-trivial randmixture integrates to 1 over the unit square", {
  skip_if_not_installed("rvinecopulib")
  vm <- randmixture(
    rvinecopulib::bicop_dist("gaussian", parameters = 0.7),
    rvinecopulib::bicop_dist("gaussian", parameters = -0.7),
    selector = function(v1, v2) pmax(v1, v2) > 0.6
  )
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3),
    randomizermod = vm
  )
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

test_that("dbsicopula()/plot() work for a randmixture built from degenerate (comonotonic/countermonotonic) copulas", {
  # regression test: this exact construction (rho = +-1) used to error inside
  # plot()/dbsicopula(), which only handled randsdvine randomizer models
  skip_if_not_installed("rvinecopulib")
  vm <- randmixture(
    rvinecopulib::bicop_dist("gaussian", parameters = 1),
    rvinecopulib::bicop_dist("gaussian", parameters = -1),
    selector = function(v1, v2) pmax(v1, v2) > 0.6
  )
  bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", parameters = 0.85), vsymmetric(), vsymmetric(),
    randomizermod = vm
  )
  d <- dbsicopula(c(0.1, 0.3, 0.5, 0.7, 0.9), c(0.2, 0.4, 0.5, 0.6, 0.8), bc)
  expect_true(all(is.finite(d) & d >= 0))

  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  expect_no_error(plot(bc))
})

test_that("dbsicopula() gives a finite, non-negative value when u1 lands exactly on a breakpoint", {
  skip_if_not_installed("rvinecopulib")
  sdv <- randsdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)))
  bc <- bsicopula(rvinecopulib::bicop_dist("gumbel", 0, 1.6), udpcosine(2), udpcosine(3),
    randomizermod = sdv
  )
  val <- dbsicopula(1 / 2, 0.4, bc)
  expect_true(is.finite(val) && val >= 0)
})

## plot() ---------------------------------------------------------------

test_that("plot() runs for both type = \"contour\" (default) and \"persp\", with and without a randsdvine", {
  skip_if_not_installed("rvinecopulib")
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })

  sdv <- randsdvine(rvinecopulib::bicop_dist("t", 0, c(0.4, 5)), rvinecopulib::bicop_dist("clayton", 0, 1.3))
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
  expect_length(captured, 8L) # median and up: 0.5, ..., 0.9, 0.95, 0.99, 0.999
  expect_true(!is.unsorted(captured, strictly = TRUE))

  plot(bc, n = 20, levels = c(0.5, 1, 2))
  expect_identical(captured, c(0.5, 1, 2))
})

test_that("plot(type = \"contour\")'s 'probs' argument controls the default levels", {
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
  plot(bc, n = 20, probs = c(0.25, 0.5, 0.75))
  expect_length(captured, 3L)
})

test_that("plot(type = \"contour\") fills by the rank-transformed density, not the raw values", {
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
    image = function(x, y, z, ...) {
      captured <<- z
      graphics::plot.new()
      graphics::plot.window(range(x), range(y))
      invisible(NULL)
    },
    .package = "udp"
  )
  n <- 20
  plot(bc, n = n)
  expect_equal(dim(captured), c(n, n))
  expect_true(all(captured > 0 & captured <= 1))
  # rank(), not the raw density, so the fill is spread across ~the full [0, 1]
  # range regardless of the density's own skew (ties -- e.g. from udpcosine's
  # structural symmetry -- are fine; rank() averages them)
  expect_equal(max(captured), 1)
  # values are floored to the lowest contour level before ranking (see next
  # test), so the minimum captured value reflects that tie, not raw density
  # spread near 0
  expect_gt(min(captured), 0.1)
})

test_that("plot(type = \"contour\") floors density below the lowest level before ranking, so the background ties to one colour", {
  skip_if_not_installed("rvinecopulib")
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  # a near-degenerate randmixture: comonotonic/countermonotonic components,
  # like the "regime-switching" vignette example, which leaves a speckle of
  # tiny but numerically nonzero density values in what should be a
  # perfectly uniform "background" region -- exactly what the floor should
  # absorb
  vm <- randmixture(
    rvinecopulib::bicop_dist("gaussian", parameters = 1),
    rvinecopulib::bicop_dist("gaussian", parameters = -1),
    selector = function(v1, v2) pmax(v1, v2) > 0.6
  )
  bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", parameters = 0.85), vsymmetric(), vsymmetric(),
    randomizermod = vm
  )

  captured <- NULL
  local_mocked_bindings(
    image = function(x, y, z, ...) {
      captured <<- z
      graphics::plot.new()
      graphics::plot.window(range(x), range(y))
      invisible(NULL)
    },
    .package = "udp"
  )
  n <- 60
  plot(bc, n = n)

  grid <- (seq_len(n) - 0.5) / n
  U <- expand.grid(u1 = grid, u2 = grid)
  dens <- matrix(dbsicopula(U$u1, U$u2, bc), n, n)
  levels <- sort(unique(stats::quantile(dens, c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999))))

  below <- dens < min(levels)
  expect_gt(sum(below), 0) # sanity check the test setup actually has a "background"
  expect_equal(length(unique(as.vector(captured[below]))), 1L)
})
