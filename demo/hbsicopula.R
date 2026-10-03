# Demonstration of pbsicopula() and hbsicopula(): the CDF, the h-functions and
# their inverses of a bivariate stochastic inversion copula.
#
# Run with  demo("hbsicopula", package = "udp")  or  source() this file.
# Takes about a minute and a half. Three base copulas (Gaussian, Gumbel, which has upper
# tail dependence, and Clayton), linear and non-linear udps, and independent
# and non-independent (randsdvine) randomizers. Everything is checked against
# simulation from rbsicopula(), which shares no code with the functions tested.

library(udp)
bd <- rvinecopulib::bicop_dist
set.seed(2026)

# ---- the models ------------------------------------------------------------

base_copulas <- list(
  gaussian = bd("gaussian", 0, 0.5),
  gumbel   = bd("gumbel", 0, 2), # upper tail dependence
  clayton  = bd("clayton", 90, 2) # rotated, so not exchangeable
)

# piecewise-linear udps: closed-form CDF and h with independent randomizers
linear_udps <- list(udp1 = vlinear(0.4), udp2 = udpzigzag(widths = c(3, 4, 3)))
# non-linear udps: tanh-sinh quadrature
nonlinear_udps <- list(udp1 = v2p(0.4, 1.5), udp2 = v3p(0.6, 1.3, 1.2))

randomizers <- list(
  independent = NULL,
  randsdvine  = randsdvine(bd("gaussian", 0, 0.7), bd("clayton", 90, 1.5), bd("gumbel", 270, 1.8))
)

models <- list()
for (b in names(base_copulas)) {
  for (u in c("linear", "nonlinear")) {
    for (r in names(randomizers)) {
      udps <- if (u == "linear") linear_udps else nonlinear_udps
      models[[paste(b, u, r, sep = " | ")]] <- bsicopula(
        base_copulas[[b]], udps$udp1, udps$udp2, randomizers[[r]]
      )
    }
  }
}

# ---- 1. CDF against simulation, all 12 models -------------------------------
# nodes = 31 keeps the double integral of the randomizer models quick; the
# independent-randomizer linear models are exact and ignore it.

cat("\n== 1. pbsicopula() against the empirical CDF of 2e5 simulated points ==\n")
pts <- cbind(c(0.2, 0.5, 0.8, 0.4), c(0.3, 0.5, 0.6, 0.9))
res <- t(vapply(names(models), function(nm) {
  m <- models[[nm]]
  sim <- rbsicopula(2e5, m)
  emp <- vapply(1:4, function(i) mean(sim[, 1] <= pts[i, 1] & sim[, 2] <= pts[i, 2]), 0)
  tm <- system.time(cdf <- pbsicopula(pts[, 1], pts[, 2], m, nodes = 31))[["elapsed"]]
  c(max_abs_diff = max(abs(cdf - emp)), seconds = tm)
}, numeric(2)))
print(round(res, 4))
cat("(Monte Carlo error of the empirical CDF is about 1e-3.)\n")

# ---- 2. h-functions: linear vs non-linear, independent vs randomizer --------
# The Gumbel base. Each panel shows h(u2 | u1) = P(U2 <= u2 | U1 = u1) for
# three values of u1 (lines) over the empirical conditional CDF of simulated
# points with U1 within 0.01 of u1 (dots).

cat("\n== 2. h-functions (plot) ==\n")
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
u2 <- seq(0, 1, length.out = 201)
for (nm in grep("^gumbel", names(models), value = TRUE)) {
  m <- models[[nm]]
  sim <- rbsicopula(4e5, m)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "u2", ylab = "h(u2 | u1)", main = sub("gumbel \\| ", "", nm), cex.main = 0.9)
  abline(0, 1, col = "grey80")
  for (k in seq_along(u1s <- c(0.1, 0.5, 0.9))) {
    lines(u2, hbsicopula(u1s[k], u2, m, nodes = 41), col = k + 1, lwd = 2)
    s <- sim[abs(sim[, 1] - u1s[k]) < 0.01, 2]
    points(sort(s)[round(seq(1, length(s), length.out = 40))], seq(0, 1, length.out = 40), col = k + 1, pch = 20, cex = 0.6)
  }
  legend("bottomright", paste("u1 =", u1s), col = 2:4, lwd = 2, bty = "n", cex = 0.8)
}
par(op)

# ---- 3. inverse h: conditional simulation -----------------------------------
# Drawing p ~ U(0,1) and applying the inverse h-function gives draws of
# U2 | U1 = u1. Compare with the points of an unconditional sample that fall
# near u1 (two-sample Kolmogorov-Smirnov test; a large p-value is agreement).

cat("\n== 3. inverse h as a conditional sampler: KS test against simulation ==\n")
u1 <- 0.3
ks <- vapply(names(models), function(nm) {
  m <- models[[nm]]
  sim <- rbsicopula(4e5, m)
  ref <- sim[abs(sim[, 1] - u1) < 0.005, 2]
  draw <- hbsicopula(rep(u1, 1500), runif(1500), m, inverse = TRUE, nodes = 41)
  suppressWarnings(ks.test(draw, ref)$p.value)
}, 0)
print(round(ks, 3))

# ---- 4. consistency: h is the derivative of the CDF -------------------------

cat("\n== 4. h(u2 | u1) against a finite difference of pbsicopula() in u1 ==\n")
e <- 1e-4
chk <- vapply(names(models), function(nm) {
  m <- models[[nm]]
  d <- (pbsicopula(0.55 + e, 0.6, m, nodes = 41) - pbsicopula(0.55 - e, 0.6, m, nodes = 41)) / (2 * e)
  c(finite_diff = d, h = hbsicopula(0.55, 0.6, m, nodes = 41))
}, numeric(2))
print(round(t(chk), 4))

# ---- 5. the CDF surface ------------------------------------------------------
# C(u1, u2) - u1 * u2 shows the departure from independence; the non-linear
# randsdvine model with the Gumbel base copula.

cat("\n== 5. C(u1, u2) - u1 u2 (plot) ==\n")
g <- (1:30 - 0.5) / 30
G <- expand.grid(u1 = g, u2 = g)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (nm in c("gumbel | linear | independent", "gumbel | nonlinear | randsdvine")) {
  z <- matrix(pbsicopula(G$u1, G$u2, models[[nm]], nodes = 25) - G$u1 * G$u2, 30)
  contour(g, g, z, nlevels = 12, xlab = "u1", ylab = "u2", main = nm, cex.main = 0.9)
}
par(op)
