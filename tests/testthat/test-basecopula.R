# The base-copula interface (R/basecopula.R): a new family plugs into
# bsicopula() by defining a class and methods for the basecopula_*()
# generics. A toy family that wraps a bicop_dist and delegates every method
# must therefore give exactly the same results as the bicop_dist itself.

bd <- function(...) rvinecopulib::bicop_dist(...)

setClass("toycop", slots = list(cop = "ANY"))
setMethod("basecopula_supported", "toycop", function(basecopula) TRUE)
setMethod("basecopula_sample", "toycop", function(basecopula, n) basecopula_sample(basecopula@cop, n))
setMethod("basecopula_density", "toycop", function(basecopula, v1, v2) basecopula_density(basecopula@cop, v1, v2))
setMethod("basecopula_cdf", "toycop", function(basecopula, v1, v2) basecopula_cdf(basecopula@cop, v1, v2))
setMethod("basecopula_h", "toycop", function(basecopula, given, x, given_var) basecopula_h(basecopula@cop, given, x, given_var))
setMethod("basecopula_hinv", "toycop", function(basecopula, given, q, given_var) basecopula_hinv(basecopula@cop, given, q, given_var))
setMethod("basecopula_fast_hinv", "toycop", function(basecopula) TRUE)

test_that("a plug-in base copula behaves like the copula it wraps", {
  skip_if_not_installed("rvinecopulib")
  real <- bd("clayton", 90, 2)
  toy <- new("toycop", cop = real)
  mk <- function(cop, rm = NULL) bsicopula(cop, vlinear(0.4), v2p(0.3, 1.4), rm)
  set.seed(1)
  u1 <- runif(20)
  u2 <- runif(20)
  a <- mk(real)
  b <- mk(toy)
  expect_equal(dbsicopula(u1, u2, b), dbsicopula(u1, u2, a))
  expect_equal(pbsicopula(u1, u2, b), pbsicopula(u1, u2, a))
  for (cv in 1:2) {
    h <- hbsicopula(u1, u2, b, cond_var = cv)
    expect_equal(h, hbsicopula(u1, u2, a, cond_var = cv))
    expect_equal(hbsicopula(u1, h, b, cond_var = cv, inverse = TRUE), hbsicopula(u1, h, a, cond_var = cv, inverse = TRUE))
  }
  set.seed(2)
  s1 <- rbsicopula(50, a)
  set.seed(2)
  expect_equal(rbsicopula(50, b), s1)
  expect_equal(pcells(b, 0.3, 0.6), pcells(a, 0.3, 0.6))
  expect_equal(pcells(b, v1 = 0.3, ngrid = 50), pcells(a, v1 = 0.3, ngrid = 50))
  # a plug-in may also be used in a randmixture
  expect_s4_class(randmixture(toy, real, function(v1, v2) v1 < v2), "randmixture")
})

test_that("the interface accepts exactly the supported base copulas", {
  expect_true(basecopula_supported(bd("gaussian", 0, 0.5)))
  expect_false(basecopula_supported("a copula"))
  expect_error(bsicopula("a copula", vlinear(0.4), vlinear(0.6)), "basecopula")
  expect_error(basecopula_density("a copula", 0.5, 0.5), "not a supported base copula")
  skip_if_not_installed("copula")
  expect_true(basecopula_supported(copula::claytonCopula(2)))
})
