# fitbsicopula(): maximum likelihood estimation of bsicopula objects. Needs
# the Suggested package rvinecopulib throughout.

pobs <- function(X) apply(X, 2, rank) / (nrow(X) + 1)

test_that("to_free() and from_free() are inverse on every kind of bound", {
  lower <- c(-1, 0, 1, -Inf)
  upper <- c(1, 1, Inf, Inf)
  x <- c(0.3, 0.9, 2.5, -4)
  expect_equal(from_free(to_free(x, lower, upper), lower, upper), x)
  # extreme unconstrained values stay strictly inside the bounds
  y <- from_free(c(-800, 800, -800, 5), lower, upper)
  expect_true(all(y[1:3] > lower[1:3]) && y[1] < upper[1] && y[2] < upper[2])
})

test_that("nudge_start() moves starting values off their bounds", {
  s <- nudge_start(c(1, 0, 0), c(1, 0, 0), c(50, 1, Inf))
  expect_true(s[1] > 1 && s[2] > 0 && s[3] > 0)
  expect_equal(nudge_start(0.3, -1, 1), 0.3)
})

test_that("basecopula_needs_vfloor() clamps every copula unbounded at some corner", {
  skip_if_not_installed("rvinecopulib")
  bd <- rvinecopulib::bicop_dist
  expect_true(basecopula_needs_vfloor(bd("gaussian", 0, 0.5)))
  expect_true(basecopula_needs_vfloor(bd("gaussian", 0, -0.5)))
  expect_true(basecopula_needs_vfloor(bd("t", 0, c(-0.3, 5))))
  expect_true(basecopula_needs_vfloor(bd("clayton", 0, 2)))
  expect_true(basecopula_needs_vfloor(bd("clayton", 90, 2)))
  expect_true(basecopula_needs_vfloor(bd("clayton", 180, 2)))
  expect_true(basecopula_needs_vfloor(bd("gumbel", 0, 2)))
  expect_true(basecopula_needs_vfloor(bd("joe", 0, 2)))
  expect_true(basecopula_needs_vfloor(bd("bb6", 0, c(2, 2))))
  expect_false(basecopula_needs_vfloor(bd("frank", 0, 5)))
  expect_false(basecopula_needs_vfloor(bd("bb8", 180, c(2, 0.7))))
  expect_false(basecopula_needs_vfloor(bd()))
})

test_that("clamp_v() clamps both ends, but only for observations inside [f, 1 - f]", {
  V <- c(0.001, 0.999, 0.5, 0.001, 0.999)
  u <- c(0.3, 0.6, 0.5, 0.005, 0.998)
  expect_equal(clamp_v(V, u, 0.01), c(0.01, 0.99, 0.5, 0.001, 0.999))
})

test_that("the clamp enters only the base copula density, and dbsicopula() is unchanged", {
  skip_if_not_installed("rvinecopulib")
  bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), vlinear(0.4), vlinear(0.6))
  u1 <- c(0.4, 0.2, 0.7)
  u2 <- c(0.6, 0.5, 0.1)
  expect_identical(dbsicopula(u1, u2, bc), dbsicopula_eval(u1, u2, bc))
  clamped <- dbsicopula_eval(u1, u2, bc, vfloor = 0.01)
  expect_equal(clamped[2:3], dbsicopula(u1[2:3], u2[2:3], bc))
  expect_equal(clamped[1], rvinecopulib::dbicop(c(0.01, 0.01), bc@basecopula))
})

test_that("fitbsicopula() validates its arguments", {
  skip_if_not_installed("rvinecopulib")
  bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), vlinear(0.4), vlinear(0.6))
  set.seed(1)
  U <- pobs(rbsicopula(100, bc))
  expect_error(fitbsicopula(U, "x"), "class 'bsicopula'")
  expect_error(fitbsicopula(cbind(U, U[, 1]), bc), "two-column")
  expect_error(fitbsicopula(rbind(U, c(0, 0.5)), bc), "boundaryadjust")
  expect_error(fitbsicopula(U, bc, pseudo = "yes"), "'pseudo' must be TRUE or FALSE")
  expect_error(fitbsicopula(U, bc, udpfix = NA), "'udpfix' must be TRUE or FALSE")
  expect_error(fitbsicopula(U, bc, se = "sandwich"), "'se' must be")
  expect_error(fitbsicopula(U, bc, se = "bootstrap", B = 1), "'B' must be")
  expect_error(fitbsicopula(U, bc, vfloor = "yes"), "'vfloor' must be")

  rmix <- randmixture(rvinecopulib::bicop_dist(), rvinecopulib::bicop_dist(),
    selector = function(v1, v2) v1 > v2
  )
  expect_error(fitbsicopula(U, bsicopula(bc@basecopula, vlinear(0.4), vlinear(0.6), rmix)),
    "randomizermod = NULL or a 'randsdvine'"
  )
  # every udp class in the package is supported; a new one without a
  # udp_fitpars() method gets the default, which points to udpfix = TRUE
  setClass("udpnewclass", contains = "udp", where = environment())
  setMethod("udptrans", "udpnewclass", function(x, u) u, where = environment())
  newudp <- new("udpnewclass")
  expect_error(
    fitbsicopula(U, bsicopula(bc@basecopula, newudp, vlinear(0.6))),
    "estimating the parameters of a 'udpnewclass' object is not yet supported; use udpfix = TRUE"
  )
  expect_s4_class(fitbsicopula(U, bsicopula(bc@basecopula, newudp, vlinear(0.6)), udpfix = TRUE),
    "fitbsicopula"
  )
  expect_error(
    fitbsicopula(U, bsicopula(rvinecopulib::bicop_dist("tll"), vlinear(0.4), vlinear(0.6))),
    "not supported for the 'tll' copula family"
  )

  skip_if_not_installed("copula")
  expect_error(
    fitbsicopula(U, bsicopula(copula::normalCopula(0.5), vlinear(0.4), vlinear(0.6))),
    "bicop_dist"
  )
})

test_that("fitbsicopula() recovers the parameters of a v-transform model", {
  skip_if_not_installed("rvinecopulib")
  truth <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.6), vlinear(0.4), v2p(0.6, 1.5))
  start <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.3), vlinear(0.5), v2p(0.5, 1))
  set.seed(1)
  U <- pobs(rbsicopula(2000, truth))
  fit <- fitbsicopula(U, start)

  expect_s4_class(fit, "fitbsicopula")
  expect_identical(
    names(coef(fit)),
    c("basecopula.rho", "udp1.delta", "udp2.delta", "udp2.kappa")
  )
  expect_equal(unname(coef(fit)), c(0.6, 0.4, 0.6, 1.5), tolerance = 0.1)
  expect_identical(fit@convergence, 0L)
  expect_identical(fit@npar, 4L)
  expect_identical(fit@nobs, 2000L)
  expect_equal(fit@vfloor, 1 / 4000)
  expect_identical(fit@se_method, "none")
  expect_true(all(is.na(fit@se)))
  expect_error(vcov(fit), "no covariance matrix")

  # the fitted object carries the estimates, on the natural scale
  expect_equal(unname(fit@bsicopula@udp2@pars), unname(coef(fit)[3:4]))
  expect_equal(fit@loglik, sum(log(dbsicopula(U[, 1], U[, 2], fit@bsicopula))), tolerance = 1e-6)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(AIC(fit), -2 * fit@loglik + 8)
  expect_output(show(fit), "udp1.delta")
})

test_that("parameter-free udps and udpfix = TRUE leave only the copula to estimate", {
  skip_if_not_installed("rvinecopulib")
  truth <- bsicopula(rvinecopulib::bicop_dist("joe", 0, 2), vsymmetric(), udpcosine(2))
  set.seed(2)
  U <- pobs(rbsicopula(500, truth))
  fit <- fitbsicopula(U, truth)
  expect_identical(names(coef(fit)), "basecopula.theta")
  expect_equal(fit@vfloor, 1 / 1000)

  bc <- bsicopula(rvinecopulib::bicop_dist("frank", 0, 3), v2p(0.45, 1.2), udpzigzag(widths = c(1, 2)))
  fit2 <- fitbsicopula(U, bc, udpfix = TRUE)
  expect_identical(names(coef(fit2)), "basecopula.theta")
  expect_identical(fit2@bsicopula@udp1@pars, bc@udp1@pars)
  expect_true(is.na(fit2@vfloor))
})

test_that("randsdvine fits estimate every parametric randomizer copula, and honour twostage", {
  skip_if_not_installed("rvinecopulib")
  bd <- rvinecopulib::bicop_dist
  truth <- bsicopula(bd("gaussian", 0, 0.5), vlinear(0.4), vlinear(0.6), randsdvine(bd("gaussian", 0, 0.7)))
  set.seed(3)
  U <- pobs(rbsicopula(1000, truth))

  fit <- fitbsicopula(U, truth)
  expect_identical(
    names(coef(fit)),
    c("basecopula.rho", "udp1.delta", "udp2.delta", "copZ1Z2_V1V2.rho")
  )
  expect_equal(unname(coef(fit)["copZ1Z2_V1V2.rho"]), 0.7, tolerance = 0.1)
  expect_s4_class(fit@stage1, "fitbsicopula")
  expect_null(fit@stage1@bsicopula@randomizermod)
  expect_null(fitbsicopula(U, truth, twostage = FALSE)@stage1)

  # the tree-2 copulas are estimated too once they have parameters to
  # estimate (the independence copula has none)
  start <- truth
  start@randomizermod <- randsdvine(bd("gaussian", 0, 0.7), bd("gaussian", 0, 0.1), bd("gaussian", 0, 0.1))
  fit3 <- fitbsicopula(U, start, udpfix = TRUE)
  expect_identical(
    names(coef(fit3)),
    c("basecopula.rho", "copZ1Z2_V1V2.rho", "copZ1V2_V1.rho", "copV1Z2_V2.rho")
  )
  expect_identical(names(coef(fitbsicopula(U, truth, udpfix = TRUE))),
    c("basecopula.rho", "copZ1Z2_V1V2.rho")
  )
})

test_that("standard errors: Hessian and parametric bootstrap", {
  skip_if_not_installed("rvinecopulib")
  truth <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.6), vlinear(0.4), vlinear(0.6))
  set.seed(4)
  U <- pobs(rbsicopula(300, truth))

  fh <- fitbsicopula(U, truth, se = "hessian")
  expect_identical(fh@se_method, "hessian")
  expect_true(all(is.finite(fh@se) & fh@se > 0))
  expect_equal(dim(vcov(fh)), c(3L, 3L))
  expect_identical(nrow(fh@boot), 0L)

  set.seed(5)
  fb <- fitbsicopula(U, truth, se = "bootstrap", B = 5)
  expect_identical(fb@se_method, "bootstrap")
  expect_identical(dim(fb@boot), c(5L, 3L))
  expect_identical(colnames(fb@boot), names(coef(fb)))
  expect_equal(fb@se, sqrt(diag(stats::cov(fb@boot))))
  expect_equal(vcov(fb), stats::cov(fb@boot))
  expect_output(show(fb), "parametric bootstrap standard errors, 5 of 5 refits")

  # the estimates don't depend on the standard-error method, and TRUE means bootstrap
  expect_equal(coef(fb), coef(fh))
  set.seed(5)
  expect_equal(fitbsicopula(U, truth, se = TRUE, B = 5)@se, fb@se)
})

test_that("zigzag breakpoints are estimated, with the number of pieces and up held fixed", {
  skip_if_not_installed("rvinecopulib")
  truth <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.6), udpid(),
    udpzigzag(widths = c(0.3, 0.45, 0.25), up = FALSE)
  )
  start <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.3), udpid(),
    udpzigzag(widths = c(1, 1, 1), up = FALSE)
  )
  set.seed(6)
  U <- pobs(rbsicopula(2000, truth))
  fit <- fitbsicopula(U, start)
  expect_identical(names(coef(fit)), c("basecopula.rho", "udp2.break1", "udp2.break2"))
  expect_equal(unname(coef(fit)), c(0.6, 0.3, 0.75), tolerance = 0.1)
  expect_false(fit@bsicopula@udp2@up)
  expect_equal(fit@bsicopula@udp2@breaks, c(0, unname(coef(fit)[2:3]), 1))

  # a one-piece zigzag has no breakpoints to estimate
  one <- bsicopula(rvinecopulib::bicop_dist("frank", 0, 3), udpid(), udpzigzag(breaks = numeric(0)))
  expect_identical(names(coef(fitbsicopula(U, one))), "basecopula.theta")
})

test_that("the zigzag's log-ratio maps are inverse and keep the breakpoints ordered", {
  maps <- udp_fitpars(udpzigzag(widths = c(1, 2, 3, 4)))$maps
  b <- c(0.1, 0.3, 0.6)
  expect_equal(maps$from_free(maps$to_free(b)), b)
  bb <- maps$from_free(c(40, -40, 3))
  expect_true(all(diff(c(0, bb, 1)) > 0))
  expect_true(all(diff(c(0, maps$nudge(c(0.001, 0.002, 0.5)), 1)) >= 0.01 - 1e-12))
})

test_that("with pseudo = FALSE extreme observations are allowed and the bootstrap skips ranking", {
  skip_if_not_installed("rvinecopulib")
  truth <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.6), vlinear(0.4), vlinear(0.6))
  set.seed(7)
  U <- rbsicopula(300, truth)
  U[1, ] <- c(1e-5, 1 - 1e-5)
  fit <- fitbsicopula(U, truth, pseudo = FALSE)
  expect_true(is.finite(fit@loglik))

  # identical seeds: the only difference is whether the bootstrap ranks
  set.seed(8)
  fp <- fitbsicopula(U, truth, pseudo = FALSE, se = "bootstrap", B = 4)
  set.seed(8)
  fq <- fitbsicopula(U, truth, pseudo = TRUE, se = "bootstrap", B = 4)
  expect_equal(coef(fp), coef(fq))
  expect_false(isTRUE(all.equal(fp@boot, fq@boot)))
})

test_that("the unit-weight chart normalizes, round-trips and keeps the pivot's sign", {
  fp <- unit_weight_fitpars(c(3, -4, 1))
  expect_equal(unname(fp$value), c(3, -4, 1) / sqrt(26))
  expect_identical(names(fp$value), c("coef1", "coef2", "coef3"))
  m <- fp$maps
  expect_length(m$to_free(fp$value), 2L)
  expect_equal(m$from_free(m$to_free(fp$value)), unname(fp$value))
  cf <- m$from_free(c(50, -80))
  expect_equal(sum(cf^2), 1)
  expect_true(cf[2] < 0) # the pivot (coef2, largest in magnitude) stays negative
  expect_identical(unit_weight_fitpars(2), no_fitpars)
})

test_that("vectorized crossings and sublevel measures agree with polynomial roots", {
  set.seed(10)
  for (r in 1:5) {
    x <- udplegendrebex(rnorm(4))
    y <- seq(x@lbound, x@ubound, length.out = 50)
    knots <- c(0, legendre_turnpoints(x@cfsD), 1)
    R <- poly_crossings(x@cfs, y, knots)
    for (i in seq_along(y)) {
      # legendre_realroots() merges a turning point's two roots; poly_crossings()
      # keeps one per branch
      expect_equal(unique(signif(R[i, !is.na(R[i, ])], 9)),
        signif(legendre_realroots(x@cfs, y[i]), 9), tolerance = 1e-6)
    }
    expect_equal(
      poly_sublevel_measure(x@cfs, y, knots),
      vapply(y, legendre_measure_bounded, 0, coef = x@cfs, lbound = x@lbound, ubound = x@ubound),
      tolerance = 1e-8
    )
  }
})

test_that("the weights of udplegendrebex and udpcosinebex objects are estimated at unit length", {
  skip_if_not_installed("rvinecopulib")
  nrm <- function(x) x / sqrt(sum(x^2))
  cf <- nrm(c(0.3, -0.8, 0.5))
  cf0 <- nrm(cf + c(0.15, 0.1, -0.2))
  for (ctor in list(udplegendrebex, udpcosinebex)) {
    truth <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.6), udpid(), ctor(cf))
    set.seed(1)
    U <- pobs(rbsicopula(1000, truth))
    fit <- fitbsicopula(U, bsicopula(truth@basecopula, udpid(), ctor(cf0)), se = "hessian")
    expect_identical(
      names(coef(fit)),
      c("basecopula.rho", "udp2.coef1", "udp2.coef2", "udp2.coef3")
    )
    expect_identical(fit@npar, 3L)
    expect_equal(AIC(fit), -2 * fit@loglik + 6)
    expect_equal(sum(coef(fit)[2:4]^2), 1)
    expect_equal(unname(coef(fit)[2:4]), cf, tolerance = 0.15)
    expect_equal(fit@bsicopula@udp2@coef, unname(coef(fit)[2:4]))
    expect_length(fit@se, 4L)
    expect_equal(dim(vcov(fit)), c(4L, 4L))
  }
})

test_that("dbsicopula() handles randsdvine models with an expansion udp", {
  skip_if_not_installed("rvinecopulib")
  bc <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.5), vlinear(0.4),
    udplegendrebex(c(0.3, -0.8, 0.5)), randsdvine(rvinecopulib::bicop_dist("gaussian", 0, 0.7))
  )
  set.seed(2)
  U <- pobs(rbsicopula(1000, bc))
  d <- dbsicopula(U[, 1], U[, 2], bc)
  expect_true(all(is.finite(d) & d > 0))
})
