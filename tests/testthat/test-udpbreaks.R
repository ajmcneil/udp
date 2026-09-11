# udpbreaks(): the internal (unexported) partition on which a udp
# transformation is piecewise continuously differentiable.

test_that("shuffle: udpbreaks() is the strip boundaries", {
  expect_equal(udpbreaks(shuffle(c(3, 1, 2))), c(0, 1 / 3, 2 / 3, 1))
  expect_equal(udpbreaks(shuffle(1)), c(0, 1))
})

test_that("udpcosine: udpbreaks() is the kink points", {
  expect_equal(udpbreaks(udpcosine(4)), (0:4) / 4)
})

test_that("vtransform: udpbreaks() is 0, delta, 1", {
  expect_equal(udpbreaks(v2p(delta = 0.37, kappa = 1.1)), c(0, 0.37, 1))
  # vsymmetric() does not store delta in `pars`; falls back to its implicit 0.5
  expect_equal(udpbreaks(vsymmetric()), c(0, 0.5, 1))
})

test_that("the default udp method assumes no break points", {
  setClass("udp_dummy_no_breaks", contains = "udp")
  expect_equal(udpbreaks(new("udp_dummy_no_breaks")), c(0, 1))
  removeClass("udp_dummy_no_breaks")
})

test_that("udplegendre: piece counts match the known values for degree 2:9", {
  # 2, 5, 6, 13, 12, 25, 20, 41 pieces respectively; degree 7 and 9 are a
  # regression check for the double-root de-duplication (see udpbreaks())
  counts <- vapply(2:9, function(d) length(udpbreaks(udplegendre(d))) - 1L, integer(1))
  expect_equal(counts, c(2L, 5L, 6L, 13L, 12L, 25L, 20L, 41L))
})

test_that("udplegendre: udpbreaks() always starts at 0, ends at 1, strictly increasing", {
  # includes degree 7+, where L_j(u) - y has an ill-conditioned double root
  # exactly at a turning point's own critical value (see udpbreaks() comment)
  for (d in 1:12) {
    b <- udpbreaks(udplegendre(d))
    expect_equal(b[1], 0)
    expect_equal(b[length(b)], 1)
    expect_true(all(diff(b) > 1e-5))
  }
  # degree 1 (T is the identity) has no interior break points
  expect_equal(udpbreaks(udplegendre(1)), c(0, 1))
})

test_that("udplegendre: every interior break point is a genuine left/right derivative mismatch", {
  eps <- 1e-4
  for (d in 2:7) {
    x <- udplegendre(d)
    b <- udpbreaks(x)
    interior <- b[b > 1e-6 & b < 1 - 1e-6]
    left <- udpderiv(x, interior - eps)
    right <- udpderiv(x, interior + eps)
    expect_true(all(abs(left - right) > 0.5))
  }
})

test_that("udplegendre: udpbreaks() agrees with the turning points and their re-crossings", {
  x <- udplegendre(4)
  tp <- legendre_turnpoints(x@cfsD)
  b <- udpbreaks(x)
  expect_true(all(vapply(tp, function(t) any(abs(b - t) < 1e-6), logical(1))))
  # degree 4 has 3 turning points but 5 interior break points: the outer two
  # turning points each have a transversal re-crossing elsewhere in [0, 1]
  expect_length(b, 7L)
  expect_length(tp, 3L)
})
