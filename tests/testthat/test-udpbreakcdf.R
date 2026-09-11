# udpbreakcdf(): P(udpsi(x, v) <= b) for each v and each break point b --
# the running sum of udpinverse()'s pre-image selection probabilities not
# exceeding b.

test_that("udpbreakcdf() matches a brute-force sum over pre-images", {
  cases <- list(
    shuffle(c(3, 1, 2), signs = c(1, -1, 1)),
    udpcosine(4),
    v2p(delta = 0.4, kappa = 1.2),
    udplegendre(4)
  )
  set.seed(1)
  v <- runif(15)
  for (x in cases) {
    b <- udpbreaks(x)
    M <- udpinverse(x, v, prob = TRUE)
    P <- attr(M, "prob")
    brute <- sapply(b, function(bi) rowSums(P * (M <= bi), na.rm = TRUE))
    expect_equal(udpbreakcdf(x, v, b), brute, ignore_attr = TRUE)
  }
})

test_that("udpbreakcdf() is 0 at the left end and 1 at the right end", {
  for (x in list(shuffle(c(2, 1, 3)), udpcosine(5),
                 v3p(delta = 0.45, kappa = 0.8, xi = 1.2), udplegendre(5))) {
    v <- seq(0.05, 0.95, length.out = 9)
    cdf <- udpbreakcdf(x, v, c(0, 1))
    expect_equal(cdf[, 1], rep(0, length(v)))
    expect_equal(cdf[, 2], rep(1, length(v)))
  }
})

test_that("udpbreakcdf() is non-decreasing across sorted break points", {
  x <- udplegendre(6)
  set.seed(2)
  v <- runif(20)
  b <- udpbreaks(x)
  cdf <- udpbreakcdf(x, v, b)
  expect_true(all(apply(cdf, 1, diff) >= -1e-9))
})

test_that("udpcosine: udpbreakcdf() at its own breaks is exactly k / degree", {
  # udpbreaks(udpcosine(d)) is exactly the piece boundaries, and every piece
  # carries equal weight 1 / d, so the CDF steps to k / d at the k-th break
  # regardless of v
  x <- udpcosine(4)
  v <- c(0.1, 0.35, 0.6, 0.9)
  b <- udpbreaks(x)
  cdf <- udpbreakcdf(x, v, b)
  expected <- matrix(rep((0:4) / 4, each = length(v)), nrow = length(v))
  expect_equal(cdf, expected, ignore_attr = TRUE)
})

test_that("vtransform: udpbreakcdf() at the fulcrum equals vdownprob()", {
  for (case in vtransform_list()) {
    v <- u_grid()
    cdf <- udpbreakcdf(case$x, v, udpbreaks(case$x))
    expect_equal(cdf[, 2], vdownprob(case$x, v), tolerance = 1e-6)
  }
})

test_that("udpbreakcdf() agrees with the empirical distribution of udpsi()", {
  skip_on_cran()
  set.seed(3)
  x <- udplegendre(4)
  v <- 0.42
  b <- udpbreaks(x)
  cdf <- udpbreakcdf(x, rep(v, 1), b)
  Z <- runif(2e4)
  draws <- udpsi(x, rep(v, length(Z)), Z)
  empirical <- vapply(b, function(bi) mean(draws <= bi), numeric(1))
  expect_equal(as.numeric(cdf), empirical, tolerance = 0.02)
})
