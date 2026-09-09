# udpinverse() enumerates every pre-image of v as an n-by-k matrix, left-packed
# and NA-padded, optionally carrying an aligned "prob" attribute that udpsi()
# samples against.

test_that("udpinverse() returns a matrix shaped by the pre-image count", {
  expect_equal(dim(udpinverse(shuffle(c(2, 1, 3)), c(0.2, 0.5, 0.9))), c(3L, 1L))
  expect_equal(dim(udpinverse(vsymmetric(), c(0.2, 0.5, 0.9))), c(3L, 2L))
  expect_equal(dim(udpinverse(udpcosine(4), c(0.2, 0.5, 0.9))), c(3L, 4L))
  expect_equal(dim(udpinverse(udpcosine(1), c(0.2, 0.5))), c(2L, 1L))
})

test_that("udpinverse() rows are the sorted pre-images, left-packed with NA", {
  x <- udpcosine(5)
  v <- c(0, 0.3, 0.7, 1)
  M <- udpinverse(x, v)
  for (i in seq_along(v)) {
    r <- M[i, ]
    filled <- r[!is.na(r)]
    expect_false(is.unsorted(filled))
    expect_true(all(diff(is.na(r)) >= 0)) # NAs only ever trail
    expect_equal(udptrans(x, filled), rep(v[i], length(filled)), tolerance = 1e-9)
  }
})

test_that("udpinverse() agrees with the class-specific primitives", {
  s <- shuffle(c(3, 1, 2), signs = c(1, -1, 1))
  v <- c(0.1, 0.5, 0.9)
  expect_equal(as.numeric(udpinverse(s, v)), shinverse(s, v))

  x <- v2p(delta = 0.4, kappa = 1.3)
  M <- udpinverse(x, v)
  expect_equal(M[, 1], vinverse(x, v))
  expect_equal(M[, 2], v + vinverse(x, v))

  cx <- udpcosine(3)
  expect_equal(
    udpinverse(cx, c(0.4, 0.6)),
    do.call(rbind, udpcosinverse(cx, c(0.4, 0.6)))
  )
})

test_that("udpinverse(prob = TRUE) attaches an aligned probability matrix", {
  objects <- list(
    shuffle(c(2, 1, 3)), vsymmetric(), v2p(delta = 0.4, kappa = 1.3),
    udpcosine(1), udpcosine(4)
  )
  v <- c(0, 0.25, 0.5, 0.75, 1)
  for (obj in objects) {
    M <- udpinverse(obj, v, prob = TRUE)
    P <- attr(M, "prob")
    expect_identical(dim(P), dim(M))
    expect_identical(is.na(P), is.na(M)) # same NA pattern
    expect_equal(rowSums(P, na.rm = TRUE), rep(1, length(v)), tolerance = 1e-12)
    expect_true(all(P[!is.na(P)] >= 0))
  }
})

test_that("udpinverse() probabilities: equal for cosine, vdownprob for vtransform", {
  P <- attr(udpinverse(udpcosine(4), c(0.3, 0.6), prob = TRUE), "prob")
  expect_equal(P, matrix(0.25, 2, 4))

  vx <- v3p(delta = 0.45, kappa = 0.8, xi = 1.2)
  v <- seq(0.1, 0.9, by = 0.1)
  P <- attr(udpinverse(vx, v, prob = TRUE), "prob")
  expect_equal(P[, 1], vdownprob(vx, v), tolerance = 1e-8)
  expect_equal(P[, 2], 1 - vdownprob(vx, v), tolerance = 1e-8)
})

test_that("udpinverse() bare call omits the prob attribute", {
  expect_null(attr(udpinverse(vsymmetric(), c(0.2, 0.6)), "prob"))
  expect_null(attr(udpinverse(udpcosine(3), c(0.2, 0.6)), "prob"))
})

test_that("udpinverse() handles empty v and drops v's attributes", {
  expect_equal(dim(udpinverse(udpcosine(3), numeric(0))), c(0L, 3L))
  expect_equal(dim(udpinverse(vsymmetric(), numeric(0))), c(0L, 2L))

  v <- ts(runif(10), frequency = 4)
  M <- udpinverse(vsymmetric(), v)
  expect_null(attr(M, "tsp"))
  expect_equal(nrow(M), 10L)
})

test_that("udpinverse() validates v for the cosine method", {
  expect_error(udpinverse(udpcosine(3), c(0.5, 1.5)), "in \\[0, 1\\]")
  expect_error(udpinverse(udpcosine(3), c(0.5, NA)), "in \\[0, 1\\]")
})

test_that("udpsi() selects branches with the udpinverse probabilities", {
  set.seed(1)
  x <- v2p(delta = 0.4, kappa = 1.6)
  v <- runif(6000)
  u <- udpsi(x, v, runif(length(v)))
  took_lower <- abs(u - vinverse(x, v)) < 1e-8
  # the fraction on the lower branch tracks the mean down-probability
  expect_equal(mean(took_lower), mean(vdownprob(x, v)), tolerance = 0.03)
})

test_that("udpsi() reproduces the equal-weight cosine selection", {
  x <- udpcosine(5)
  v <- seq(0.05, 0.95, length.out = 40)
  roots <- do.call(rbind, udpcosinverse(x, v))
  for (z in c(1e-9, 0.5, 1)) {
    j <- min(ceiling(z * 5), 5)
    expect_equal(udpsi(x, v, rep(z, length(v))), roots[, j])
  }
})
