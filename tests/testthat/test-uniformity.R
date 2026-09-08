# The defining property of the package: a v-transform maps a standard uniform
# to a standard uniform, and its stochastic inverse maps back to a standard
# uniform. These use a fixed seed so the Kolmogorov-Smirnov p-values are
# deterministic; they guard against gross violations, not tail behaviour.

test_that("v-transforms preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (case in vtransform_list()) {
    v <- vtrans(case$x, u)
    expect_gte(suppressWarnings(stats::ks.test(v, "punif")$p.value), 0.01)
  }
})

test_that("the stochastic inverse recovers a uniform input", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (case in vtransform_list()) {
    back <- vsi(case$x, vtrans(case$x, u))
    expect_gte(suppressWarnings(stats::ks.test(back, "punif")$p.value), 0.01)
  }
})

test_that("vsi() inverts vtrans() pointwise for a uniform input", {
  skip_on_cran()
  set.seed(1)
  u <- runif(5000)
  for (case in vtransform_list()) {
    back <- vsi(case$x, vtrans(case$x, u), Z = rep(0, length(u)))
    # Z = 0 always takes the lower pre-image; recovers u wherever u <= delta
    lower <- u <= case$delta
    expect_equal(back[lower], u[lower], tolerance = 1e-6)
  }
})
