# The order of the variables of every copula in a randsdvine model. The
# D-vine is (Z1, V1, V2, Z2) and each pair-copula is specified in the order its
# variables appear there: basecopula (V1, V2), copZ1V2_V1 (Z1, V2),
# copV1Z2_V2 (V1, Z2), copZ1Z2_V1V2 (Z1, Z2). Transposing a copula that is not
# exchangeable -- a 90 or 270 degree rotation, where it swaps the two -- gives
# another copula, so a swapped argument is an error that exchangeable
# copulas never reveal.
#
# Every check below compares the sampler with pbicop() alone, never with
# hbicop(), and conditional distributions are obtained by numerical
# differentiation of pbicop(), so that a mistake in the package's use of the
# h-functions cannot be repeated in the test. A test that mirrors the
# implementation shares its convention and would pass.

# Absolute, not relative, closeness: the simulated probabilities below are
# small (~0.03-0.07), and Monte Carlo error is absolute (~0.0008 at 1e5 draws).
expect_close <- function(x, y, tol) {
  expect_lt(max(abs(x - y)), tol)
}

# P(B <= b | A = a) (given = "first") or P(A <= a | B = b) (given = "second")
# for (A, B) ~ cop, by differentiating the copula CDF.
cond_cdf <- function(a, b, cop, given) {
  d <- 1e-6
  if (given == "first") {
    lo <- pmax(a - d, 0)
    hi <- pmin(a + d, 1)
    (rvinecopulib::pbicop(cbind(hi, b), cop) - rvinecopulib::pbicop(cbind(lo, b), cop)) / (hi - lo)
  } else {
    lo <- pmax(b - d, 0)
    hi <- pmin(b + d, 1)
    (rvinecopulib::pbicop(cbind(a, hi), cop) - rvinecopulib::pbicop(cbind(a, lo), cop)) / (hi - lo)
  }
}

# Simulate the randomizers for carriers drawn from the base copula, and the
# tree-1 conditional CDF levels e21 = P(V2 <= v2 | V1 = v1) and
# e12 = P(V1 <= v1 | V2 = v2) of each draw.
simulate_sdvine <- function(base, rm, N = 1e5) {
  V <- rvinecopulib::rbicop(N, base)
  Z <- randsdvine_sample(V[, 1], V[, 2], base, rm)
  list(
    Z1 = Z[, "Z1"], Z2 = Z[, "Z2"],
    e21 = cond_cdf(V[, 1], V[, 2], base, "first"),
    e12 = cond_cdf(V[, 1], V[, 2], base, "second")
  )
}

# The two readings of a copula at (a, b): in its own order and transposed.
# A test is only informative where they differ.
order_gap <- function(cop, a, b) {
  abs(rvinecopulib::pbicop(c(a, b), cop) - rvinecopulib::pbicop(c(b, a), cop))
}

test_that("every pair-copula of a randsdvine model is used in D-vine order", {
  skip_if_not_installed("rvinecopulib")
  bd <- rvinecopulib::bicop_dist
  pb <- rvinecopulib::pbicop
  tol <- 0.006
  set.seed(1)
  for (base in list(bd("clayton", 90, 2), bd("gumbel", 270, 2))) {
    for (rot in c(90, 270)) {
      cop <- bd("clayton", rot, 3)

      # tree 2, (Z1, V2) given V1: Z1 and e21 = F(V2 | V1) are both uniform
      # and their joint law is copZ1V2_V1 in that order
      rm <- randsdvine(bd("gaussian", 0, 0.3), cop, bd())
      S <- simulate_sdvine(base, rm)
      expect_gt(order_gap(cop, 0.3, 0.6), 0.02)
      expect_close(mean(S$Z1 <= 0.3 & S$e21 <= 0.6), pb(c(0.3, 0.6), cop), tol)

      # tree 2, (V1, Z2) given V2: e12 = F(V1 | V2) and Z2
      rm <- randsdvine(bd("gaussian", 0, 0.3), bd(), cop)
      S <- simulate_sdvine(base, rm)
      expect_close(mean(S$e12 <= 0.3 & S$Z2 <= 0.6), pb(c(0.3, 0.6), cop), tol)

      # tree 3, (Z1, Z2) given (V1, V2): the conditional CDFs
      # F(Z1 | V1, V2) and F(Z2 | V1, V2), obtained from the tree-2 copulas by
      # differentiation, have joint law copZ1Z2_V1V2 in that order
      rm <- randsdvine(cop, bd("clayton", 90, 2), bd("gumbel", 270, 2))
      S <- simulate_sdvine(base, rm)
      a <- cond_cdf(S$Z1, S$e21, rm@copZ1V2_V1, "second")
      b <- cond_cdf(S$e12, S$Z2, rm@copV1Z2_V2, "first")
      expect_close(mean(a <= 0.3 & b <= 0.6), pb(c(0.3, 0.6), cop), tol)
    }
  }
})

test_that("randomizer_cdf() is the conditional CDF of the sampled randomizers, in every slot", {
  skip_if_not_installed("rvinecopulib")
  bd <- rvinecopulib::bicop_dist
  set.seed(2)
  N <- 1e5
  base <- bd("clayton", 90, 2)
  rms <- list(
    randsdvine(bd("clayton", 90, 3), bd("clayton", 270, 2), bd("gumbel", 270, 2)),
    randsdvine(bd("clayton", 270, 3), bd("gumbel", 90, 2), bd("clayton", 270, 2))
  )
  z <- rbind(c(0.2, 0.3), c(0.5, 0.5), c(0.8, 0.4), c(0.6, 0.9))
  for (rm in rms) {
    for (v in list(c(0.3, 0.6), c(0.7, 0.2))) {
      Z <- randsdvine_sample(rep(v[1], N), rep(v[2], N), base, rm)
      F <- randomizer_cdf(v[1], v[2], base, rm)
      for (i in seq_len(nrow(z))) {
        expect_close(F(z[i, 1], z[i, 2]), mean(Z[, 1] <= z[i, 1] & Z[, 2] <= z[i, 2]), 0.006)
      }
    }
  }
})

test_that("the density and the sampler agree when every slot is an unsymmetric copula", {
  skip_if_not_installed("rvinecopulib")
  bd <- rvinecopulib::bicop_dist
  bc <- bsicopula(
    bd("clayton", 90, 2), vlinear(0.4), vlinear(0.6),
    randsdvine(bd("clayton", 270, 3), bd("clayton", 90, 3), bd("gumbel", 270, 2))
  )
  set.seed(3)
  U <- rbsicopula(1e5, bc)
  # box probabilities from the density by the midpoint rule on a 200 x 200
  # grid: the error is first order in the grid spacing (~0.0035 here, halving
  # as the grid is refined) because the density is steep near the corners
  g <- (seq_len(200) - 0.5) / 200
  G <- expand.grid(u1 = g, u2 = g)
  cuts <- seq(0, 1, 0.25)
  dens <- tapply(dbsicopula(G$u1, G$u2, bc), list(cut(G$u1, cuts), cut(G$u2, cuts)), mean) / 16
  samp <- table(cut(U[, 1], cuts), cut(U[, 2], cuts)) / nrow(U)
  expect_lt(max(abs(dens - samp)), 0.008)
  expect_gt(suppressWarnings(stats::ks.test(U[, 1], "punif")$p.value), 1e-3)
  expect_gt(suppressWarnings(stats::ks.test(U[, 2], "punif")$p.value), 1e-3)
})

test_that("transposing an unsymmetric slot changes the model", {
  skip_if_not_installed("rvinecopulib")
  # the tests above have power only if a transposed slot is distinguishable:
  # swapping 90 and 270 in copZ1V2_V1 must move the joint law of (Z1, V2 | V1)
  bd <- rvinecopulib::bicop_dist
  base <- bd("clayton", 90, 2)
  set.seed(4)
  S90 <- simulate_sdvine(base, randsdvine(bd("gaussian", 0, 0.3), bd("clayton", 90, 3), bd()))
  S270 <- simulate_sdvine(base, randsdvine(bd("gaussian", 0, 0.3), bd("clayton", 270, 3), bd()))
  expect_gt(abs(mean(S90$Z1 <= 0.3 & S90$e21 <= 0.6) - mean(S270$Z1 <= 0.3 & S270$e21 <= 0.6)), 0.02)
})
